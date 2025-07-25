
import numpy as np
from scipy.integrate import solve_ivp
import time

# Make sure to install scipy:
# pip install scipy

class KerrSpacetime:
    """
    A class to represent Kerr spacetime and calculate particle geodesics.
    Includes both single-particle and optimized batched integration methods.
    """
    def __init__(self, M=1.0, a=0.5):
        if not (0 <= a <= M):
            raise ValueError("Spin parameter 'a' must be between 0 and M.")
        self.M = M
        self.a = a

    # --- Core Component Functions (now vectorized to handle batches) ---

    def _metric_components(self, r, theta):
        """
        Calculates covariant metric components g_μν for a batch of points.

        Args:
            r (np.ndarray): Array of radial positions, shape (N,).
            theta (np.ndarray): Array of polar angles, shape (N,).

        Returns:
            np.ndarray: Batch of metric tensors, shape (N, 4, 4).
        """
        M, a = self.M, self.a
        sin2 = np.sin(theta)**2
        cos2 = np.cos(theta)**2

        Sigma = r**2 + a**2 * cos2
        Delta = r**2 - 2*M*r + a**2

        # Handle division by zero safely
        g = np.zeros((len(r), 4, 4))
        valid_mask = (Sigma > 1e-10) & (np.abs(Delta) > 1e-10)

        s = Sigma[valid_mask]
        d = Delta[valid_mask]

        g[valid_mask, 0, 0] = -(1 - 2*M*r[valid_mask] / s)
        g[valid_mask, 0, 3] = -(2*M*a*r[valid_mask]*sin2[valid_mask] / s)
        g[valid_mask, 3, 0] = g[valid_mask, 0, 3] # Symmetric
        g[valid_mask, 1, 1] = s / d
        g[valid_mask, 2, 2] = s
        g[valid_mask, 3, 3] = (r[valid_mask]**2 + a**2 + 2*M*a**2*r[valid_mask]*sin2[valid_mask] / s) * sin2[valid_mask]

        return g

    def _inverse_metric_components(self, r, theta):
        """Calculates contravariant metric g^μν for a batch of points."""
        g_batch = self._metric_components(r, theta)
        # np.linalg.inv can operate on a stack of matrices
        return np.linalg.inv(g_batch)

    def _batched_christoffel_symbols(self, r, theta):
        """Calculates Christoffel symbols Γ^μ_αβ for a batch of points."""
        # Use numerical differentiation on the batched metric function
        dr = 1e-6
        d_theta = 1e-6

        g_inv = self._inverse_metric_components(r, theta) # Shape (N, 4, 4)

        g_point = self._metric_components(r, theta)
        g_r_plus = self._metric_components(r + dr, theta)
        g_theta_plus = self._metric_components(r, theta + d_theta)

        dg_dr = (g_r_plus - g_point) / dr # Shape (N, 4, 4)
        dg_d_theta = (g_theta_plus - g_point) / d_theta # Shape (N, 4, 4)

        # dg[n, c, a, b] = ∂g_ab / ∂x^c for particle n
        # Transpose to get the correct shape: (N, 4, 4, 4)
        dg = np.stack([
            np.zeros_like(g_point), # dg/dt
            dg_dr,                  # dg/dr
            dg_d_theta,             # dg/dtheta
            np.zeros_like(g_point)  # dg/dphi
        ], axis=1)  # Shape (N, 4, 4, 4)

        # Use einsum for efficient tensor contraction
        # Γ^k_ij = 0.5 * g^kl * (∂g_li/∂x^j + ∂g_lj/∂x^i - ∂g_ij/∂x^l)
        # n=batch, k,l,i,j=indices
        Gamma = 0.5 * np.einsum('nkl,njli->nkij', g_inv, dg)
        Gamma += 0.5 * np.einsum('nkl,nilj->nkij', g_inv, dg)
        Gamma -= 0.5 * np.einsum('nkl,nlij->nkij', g_inv, dg)

        return Gamma # Shape (N, 4, 4, 4)

    def batched_geodesic_equation(self, tau, flat_state):
        """
        The vectorized geodesic equation for a batch of N particles.

        Args:
            tau: Proper time (unused, required by solver).
            flat_state (array): A flattened 1D array of shape (8*N,)
                                representing the states of all particles.

        Returns:
            array: A flattened 1D array of derivatives, shape (8*N,).
        """
        num_particles = len(flat_state) // 8
        state = flat_state.reshape((num_particles, 8))

        positions = state[:, :4] # Shape (N, 4)
        velocities = state[:, 4:] # Shape (N, 4)

        r, theta = positions[:, 1], positions[:, 2]

        # Calculate Christoffel symbols for the entire batch
        Gamma_batch = self._batched_christoffel_symbols(r, theta) # Shape (N, 4, 4, 4)

        # Calculate the 4-acceleration for the entire batch using einsum
        # a^i = -Γ^i_jk u^j u^k
        # n=batch, i,j,k=indices
        accelerations = -np.einsum('nijk,nj,nk->ni', Gamma_batch, velocities, velocities)

        # The derivative is [velocities, accelerations]
        derivatives = np.concatenate([velocities, accelerations], axis=1)

        # Return the flattened array for the solver
        return derivatives.flatten()

    def integrate_batch_trajectory(self, initial_positions, initial_velocities, dtau):
        """
        Integrates a batch of particle trajectories for a step dτ.

        Args:
            initial_positions (np.ndarray): Shape (N, 4).
            initial_velocities (np.ndarray): Shape (N, 4).
            dtau (float): Proper time step.

        Returns:
            tuple: (final_positions, final_velocities) as (N, 4) arrays.
        """
        num_particles = initial_positions.shape[0]
        if num_particles == 0:
            return np.array([]), np.array([])

        # Flatten the initial states into a single 1D vector
        initial_state_flat = np.concatenate([initial_positions, initial_velocities], axis=1).flatten()

        solution = solve_ivp(
            self.batched_geodesic_equation,
            t_span=[0, dtau],
            y0=initial_state_flat,
            method='RK45',
            t_eval=[dtau],
            rtol=1e-8,
            atol=1e-10
        )

        if not solution.success:
            raise RuntimeError(f"Integration failed: {solution.message}")

        final_state_flat = solution.y[:, -1]

        # Reshape the final flat state back into a (N, 8) array
        final_state = final_state_flat.reshape((num_particles, 8))

        final_positions = final_state[:, :4]
        final_velocities = final_state[:, 4:]

        return final_positions, final_velocities

    def integrate_one_step(self, positions, velocities, dtau):
        """
        Convenience method to integrate a batch of particles for one time step.
        
        Args:
            positions (np.ndarray): Initial positions, shape (N, 4).
            velocities (np.ndarray): Initial velocities, shape (N, 4).
            dtau (float): Proper time step.
            
        Returns:
            tuple: (new_positions, new_velocities) as (N, 4) arrays.
        """
        return self.integrate_batch_trajectory(positions, velocities, dtau)

# ==============================================================================
# USAGE EXAMPLE
# ==============================================================================
if __name__ == '__main__':
    kerr_bh = KerrSpacetime(M=1.0, a=0.9) # High spin

    # --- Set up a BATCH of 16 particles in a ring in the equatorial plane ---
    num_particles = 16
    r_initial = 15.0
    theta_initial = np.pi / 2

    # Create a ring of particles at different phi values
    phi_values = np.linspace(0, 2*np.pi, num_particles, endpoint=False)

    initial_positions = np.zeros((num_particles, 4))
    initial_positions[:, 1] = r_initial
    initial_positions[:, 2] = theta_initial
    initial_positions[:, 3] = phi_values

    # --- Set up initial velocities for a stable(ish) circular orbit ---
    # This requires finding the correct orbital velocity at this radius
    M, a = kerr_bh.M, kerr_bh.a
    # For circular orbits in the equatorial plane of Kerr spacetime
    # ω = ±M^(1/2) / (r^(3/2) ± a*M^(1/2))
    # Using the prograde orbit (+ sign)
    omega = np.sqrt(M) / (r_initial**1.5 + a * np.sqrt(M))
    dphi_dt = omega
    dr_dt = 0

    # Normalize the 4-velocity for each particle (it's the same for all in this case)
    g = kerr_bh._metric_components(initial_positions[:, 1], initial_positions[:, 2]) # Shape (N, 4, 4)
    norm_factor_sq = g[:,0,0] + 2*g[:,0,3]*dphi_dt + g[:,3,3]*dphi_dt**2
    
    # Ensure we don't take the square root of a positive number (should be negative)
    if np.any(norm_factor_sq >= 0):
        print("Warning: Found timelike normalization issue, adjusting...")
        # For particles with issues, use a smaller orbital velocity
        problematic_mask = norm_factor_sq >= 0
        omega_adjusted = omega * 0.8  # Reduce by 20%
        dphi_dt_adjusted = omega_adjusted
        norm_factor_sq_adjusted = g[:,0,0] + 2*g[:,0,3]*dphi_dt_adjusted + g[:,3,3]*dphi_dt_adjusted**2
        norm_factor_sq[problematic_mask] = norm_factor_sq_adjusted[problematic_mask]
        dphi_dt = dphi_dt_adjusted
    
    dt_dtau_batch = np.sqrt(-1 / norm_factor_sq)

    initial_velocities = np.zeros((num_particles, 4))
    initial_velocities[:, 0] = dt_dtau_batch
    initial_velocities[:, 3] = dphi_dt * dt_dtau_batch

    # --- Run the Batched Integration ---
    proper_time_step = 10.0 # Let them orbit for a bit

    print(f"--- Integrating Trajectory for a BATCH of {num_particles} Particles ---")
    print(f"Starting at radius r={r_initial} in the equatorial plane.")
    print(f"Integrating for a proper time dτ = {proper_time_step}")

    start_time = time.time()
    final_positions, final_velocities = kerr_bh.integrate_batch_trajectory(
        initial_positions, initial_velocities, proper_time_step
    )
    end_time = time.time()

    print(f"\nBatched computation finished in {end_time - start_time:.4f} seconds.")

    print("\n--- Final Positions (t, r, θ, φ) ---")
    print(np.round(final_positions, 4))

    # Check if the particles have moved in phi as expected
    delta_phi = final_positions[:, 3] - initial_positions[:, 3]
    print(f"\nAverage change in φ: {np.mean(delta_phi):.4f} radians")

    # --- Test single step integration ---
    print("\n--- Testing Single Step Integration ---")
    step_size = 0.1
    pos_step, vel_step = kerr_bh.integrate_one_step(initial_positions, initial_velocities, step_size)
    print(f"After one step of dτ = {step_size}:")
    print(f"Average φ change: {np.mean(pos_step[:, 3] - initial_positions[:, 3]):.6f} radians")
    print(f"r remains constant: {np.allclose(pos_step[:, 1], initial_positions[:, 1], rtol=1e-3)}")
    print(f"θ remains constant: {np.allclose(pos_step[:, 2], initial_positions[:, 2], rtol=1e-3)}")

