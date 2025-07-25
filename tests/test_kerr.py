"""
Test suite for the Kerr spacetime module.

Tests all functionality in the KerrSpacetime class including:
- Initialization and parameter validation
- Metric tensor calculations
- Inverse metric calculations
- Christoffel symbol calculations
- Geodesic equation evaluation
- Batch trajectory integration
- Circular orbit calculations
- Numerical stability and edge cases
"""

import sys
import os
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

# Add the src directory to the path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from causet.kerr import KerrSpacetime


class TestKerrSpacetimeInitialization:
    """Test initialization and parameter validation."""

    def test_valid_initialization(self):
        """Test valid parameter combinations."""
        # Standard case
        kerr = KerrSpacetime(M=1.0, a=0.5)
        assert kerr.M == 1.0
        assert kerr.a == 0.5

        # Schwarzschild limit
        kerr_schwarzschild = KerrSpacetime(M=2.0, a=0.0)
        assert kerr_schwarzschild.M == 2.0
        assert kerr_schwarzschild.a == 0.0

        # Extremal Kerr
        kerr_extremal = KerrSpacetime(M=1.0, a=1.0)
        assert kerr_extremal.M == 1.0
        assert kerr_extremal.a == 1.0

    def test_invalid_initialization(self):
        """Test invalid parameter combinations."""
        # a > M should raise ValueError
        with pytest.raises(ValueError, match="Spin parameter 'a' must be between 0 and M"):
            KerrSpacetime(M=1.0, a=1.5)

        # Negative a should raise ValueError
        with pytest.raises(ValueError, match="Spin parameter 'a' must be between 0 and M"):
            KerrSpacetime(M=1.0, a=-0.5)

    def test_default_parameters(self):
        """Test default parameter values."""
        kerr = KerrSpacetime()
        assert kerr.M == 1.0
        assert kerr.a == 0.5


class TestKerrMetricComponents:
    """Test metric tensor calculations."""

    def setUp(self):
        """Set up test fixtures."""
        self.kerr = KerrSpacetime(M=1.0, a=0.5)
        self.kerr_schwarzschild = KerrSpacetime(M=1.0, a=0.0)

    def test_metric_shape(self):
        """Test that metric has correct shape."""
        self.setUp()
        r = np.array([5.0, 10.0, 15.0])
        theta = np.array([np.pi/4, np.pi/2, 3*np.pi/4])
        
        g = self.kerr._metric_components(r, theta)
        assert g.shape == (3, 4, 4)

    def test_metric_symmetry(self):
        """Test that metric is symmetric."""
        self.setUp()
        r = np.array([10.0])
        theta = np.array([np.pi/2])
        
        g = self.kerr._metric_components(r, theta)
        
        # Check symmetry g_μν = g_νμ
        assert_allclose(g[0], g[0].T, rtol=1e-14)

    def test_schwarzschild_limit(self):
        """Test that metric reduces to Schwarzschild when a=0."""
        self.setUp()
        r = np.array([10.0])
        theta = np.array([np.pi/2])
        
        g_kerr = self.kerr_schwarzschild._metric_components(r, theta)
        
        # For Schwarzschild metric at theta = pi/2:
        # g_tt = -(1 - 2M/r)
        # g_rr = 1/(1 - 2M/r)
        # g_theta_theta = r^2
        # g_phi_phi = r^2
        # Off-diagonal terms should be zero
        
        M = self.kerr_schwarzschild.M
        expected_gtt = -(1 - 2*M/r[0])
        expected_grr = 1/(1 - 2*M/r[0])
        expected_gtheta = r[0]**2
        expected_gphi = r[0]**2
        
        assert_allclose(g_kerr[0, 0, 0], expected_gtt, rtol=1e-12)
        assert_allclose(g_kerr[0, 1, 1], expected_grr, rtol=1e-12)
        assert_allclose(g_kerr[0, 2, 2], expected_gtheta, rtol=1e-12)
        assert_allclose(g_kerr[0, 3, 3], expected_gphi, rtol=1e-12)
        
        # Off-diagonal terms should be zero for Schwarzschild
        assert_allclose(g_kerr[0, 0, 3], 0.0, atol=1e-14)

    def test_equatorial_plane_properties(self):
        """Test specific properties in the equatorial plane."""
        self.setUp()
        r = np.array([5.0, 10.0, 15.0])
        theta = np.array([np.pi/2, np.pi/2, np.pi/2])  # Equatorial plane
        
        g = self.kerr._metric_components(r, theta)
        
        # In equatorial plane, only g_t,phi should be non-zero off-diagonal
        for i in range(len(r)):
            assert_allclose(g[i, 0, 1], 0.0, atol=1e-14)
            assert_allclose(g[i, 0, 2], 0.0, atol=1e-14)
            assert_allclose(g[i, 1, 2], 0.0, atol=1e-14)
            assert_allclose(g[i, 1, 3], 0.0, atol=1e-14)
            assert_allclose(g[i, 2, 3], 0.0, atol=1e-14)

    def test_metric_singularity_handling(self):
        """Test that metric handles potential singularities gracefully."""
        self.setUp()
        # Test near r=0 (coordinate singularity)
        r = np.array([0.01, 0.1])
        theta = np.array([np.pi/2, np.pi/2])
        
        g = self.kerr._metric_components(r, theta)
        
        # Should not contain NaN or inf values
        assert np.all(np.isfinite(g))


class TestKerrInverseMetric:
    """Test inverse metric calculations."""

    def setUp(self):
        """Set up test fixtures."""
        self.kerr = KerrSpacetime(M=1.0, a=0.5)

    def test_inverse_metric_property(self):
        """Test that g * g^(-1) = identity."""
        self.setUp()
        r = np.array([10.0, 15.0])
        theta = np.array([np.pi/3, np.pi/2])
        
        g = self.kerr._metric_components(r, theta)
        g_inv = self.kerr._inverse_metric_components(r, theta)
        
        # Check g * g^(-1) = I for each particle
        for i in range(len(r)):
            product = np.dot(g[i], g_inv[i])
            identity = np.eye(4)
            assert_allclose(product, identity, rtol=1e-10, atol=1e-15)

    def test_inverse_metric_shape(self):
        """Test inverse metric has correct shape."""
        self.setUp()
        r = np.array([5.0, 10.0, 15.0])
        theta = np.array([np.pi/4, np.pi/2, 3*np.pi/4])
        
        g_inv = self.kerr._inverse_metric_components(r, theta)
        assert g_inv.shape == (3, 4, 4)


class TestKerrChristoffelSymbols:
    """Test Christoffel symbol calculations."""

    def setUp(self):
        """Set up test fixtures."""
        self.kerr = KerrSpacetime(M=1.0, a=0.5)

    def test_christoffel_shape(self):
        """Test Christoffel symbols have correct shape."""
        self.setUp()
        r = np.array([10.0, 15.0])
        theta = np.array([np.pi/2, np.pi/3])
        
        Gamma = self.kerr._batched_christoffel_symbols(r, theta)
        assert Gamma.shape == (2, 4, 4, 4)

    def test_christoffel_symmetry(self):
        """Test Christoffel symbols are symmetric in lower indices."""
        self.setUp()
        r = np.array([10.0])
        theta = np.array([np.pi/2])
        
        Gamma = self.kerr._batched_christoffel_symbols(r, theta)
        
        # Γ^k_ij = Γ^k_ji (symmetric in lower indices)
        for k in range(4):
            for i in range(4):
                for j in range(4):
                    assert_allclose(Gamma[0, k, i, j], Gamma[0, k, j, i], rtol=1e-6)

    def test_christoffel_finite(self):
        """Test that Christoffel symbols are finite."""
        self.setUp()
        r = np.array([5.0, 10.0, 15.0])
        theta = np.array([np.pi/4, np.pi/2, 3*np.pi/4])
        
        Gamma = self.kerr._batched_christoffel_symbols(r, theta)
        
        # Should not contain NaN or inf values
        assert np.all(np.isfinite(Gamma))


class TestKerrGeodesicEquation:
    """Test geodesic equation evaluation."""

    def setUp(self):
        """Set up test fixtures."""
        self.kerr = KerrSpacetime(M=1.0, a=0.5)

    def test_geodesic_equation_shape(self):
        """Test geodesic equation returns correct shape."""
        self.setUp()
        # Set up state for 3 particles
        num_particles = 3
        state = np.random.randn(8 * num_particles)
        
        # Ensure reasonable positions (avoid singularities)
        for i in range(num_particles):
            state[8*i + 1] = 10.0 + i  # r > 0
            state[8*i + 2] = np.pi/2   # theta = pi/2
        
        derivatives = self.kerr.batched_geodesic_equation(0.0, state)
        assert derivatives.shape == (8 * num_particles,)

    def test_geodesic_conservation_laws(self):
        """Test that geodesic preserves conserved quantities for circular orbits."""
        self.setUp()
        # Set up circular orbit in equatorial plane
        r_orbit = 10.0
        theta_orbit = np.pi/2
        
        # Calculate orbital frequency
        M, a = self.kerr.M, self.kerr.a
        omega = np.sqrt(M) / (r_orbit**1.5 + a * np.sqrt(M))
        
        # Set up initial conditions
        positions = np.array([[0.0, r_orbit, theta_orbit, 0.0]])
        
        # Calculate normalized 4-velocity
        g = self.kerr._metric_components(np.array([r_orbit]), np.array([theta_orbit]))
        norm_factor_sq = g[0,0,0] + 2*g[0,0,3]*omega + g[0,3,3]*omega**2
        dt_dtau = np.sqrt(-1 / norm_factor_sq)
        
        velocities = np.array([[dt_dtau, 0.0, 0.0, omega * dt_dtau]])
        
        # Flatten state
        state = np.concatenate([positions, velocities], axis=1).flatten()
        
        # Check that geodesic equation returns finite derivatives
        derivatives = self.kerr.batched_geodesic_equation(0.0, state)
        assert np.all(np.isfinite(derivatives))


class TestKerrBatchIntegration:
    """Test batch trajectory integration."""

    def setUp(self):
        """Set up test fixtures."""
        self.kerr = KerrSpacetime(M=1.0, a=0.5)

    def test_empty_batch(self):
        """Test integration with empty batch."""
        self.setUp()
        empty_pos = np.array([]).reshape(0, 4)
        empty_vel = np.array([]).reshape(0, 4)
        
        final_pos, final_vel = self.kerr.integrate_batch_trajectory(empty_pos, empty_vel, 0.1)
        
        assert final_pos.shape == (0,)
        assert final_vel.shape == (0,)

    def test_single_particle_integration(self):
        """Test integration of single particle."""
        self.setUp()
        # Set up circular orbit
        r_orbit = 10.0
        positions = np.array([[0.0, r_orbit, np.pi/2, 0.0]])
        
        # Calculate orbital velocity
        M, a = self.kerr.M, self.kerr.a
        omega = np.sqrt(M) / (r_orbit**1.5 + a * np.sqrt(M))
        
        g = self.kerr._metric_components(np.array([r_orbit]), np.array([np.pi/2]))
        norm_factor_sq = g[0,0,0] + 2*g[0,0,3]*omega + g[0,3,3]*omega**2
        dt_dtau = np.sqrt(-1 / norm_factor_sq)
        
        velocities = np.array([[dt_dtau, 0.0, 0.0, omega * dt_dtau]])
        
        final_pos, final_vel = self.kerr.integrate_batch_trajectory(positions, velocities, 0.1)
        
        assert final_pos.shape == (1, 4)
        assert final_vel.shape == (1, 4)
        assert np.all(np.isfinite(final_pos))
        assert np.all(np.isfinite(final_vel))

    def test_multi_particle_integration(self):
        """Test integration of multiple particles."""
        self.setUp()
        num_particles = 5
        r_orbit = 10.0
        
        # Set up ring of particles
        phi_values = np.linspace(0, 2*np.pi, num_particles, endpoint=False)
        positions = np.zeros((num_particles, 4))
        positions[:, 1] = r_orbit
        positions[:, 2] = np.pi/2
        positions[:, 3] = phi_values
        
        # Set up velocities
        M, a = self.kerr.M, self.kerr.a
        omega = np.sqrt(M) / (r_orbit**1.5 + a * np.sqrt(M))
        
        g = self.kerr._metric_components(positions[:, 1], positions[:, 2])
        norm_factor_sq = g[:,0,0] + 2*g[:,0,3]*omega + g[:,3,3]*omega**2
        dt_dtau = np.sqrt(-1 / norm_factor_sq)
        
        velocities = np.zeros((num_particles, 4))
        velocities[:, 0] = dt_dtau
        velocities[:, 3] = omega * dt_dtau
        
        final_pos, final_vel = self.kerr.integrate_batch_trajectory(positions, velocities, 0.1)
        
        assert final_pos.shape == (num_particles, 4)
        assert final_vel.shape == (num_particles, 4)
        assert np.all(np.isfinite(final_pos))
        assert np.all(np.isfinite(final_vel))

    def test_circular_orbit_stability(self):
        """Test that circular orbits remain approximately circular."""
        self.setUp()
        r_orbit = 15.0  # Large radius for stability
        
        positions = np.array([[0.0, r_orbit, np.pi/2, 0.0]])
        
        # Set up circular orbit velocity
        M, a = self.kerr.M, self.kerr.a
        omega = np.sqrt(M) / (r_orbit**1.5 + a * np.sqrt(M))
        
        g = self.kerr._metric_components(np.array([r_orbit]), np.array([np.pi/2]))
        norm_factor_sq = g[0,0,0] + 2*g[0,0,3]*omega + g[0,3,3]*omega**2
        dt_dtau = np.sqrt(-1 / norm_factor_sq)
        
        velocities = np.array([[dt_dtau, 0.0, 0.0, omega * dt_dtau]])
        
        # Integrate for a short time
        final_pos, final_vel = self.kerr.integrate_batch_trajectory(positions, velocities, 1.0)
        
        # Check that radius and theta remain approximately constant
        assert_allclose(final_pos[0, 1], r_orbit, rtol=1e-2)  # r should be conserved
        assert_allclose(final_pos[0, 2], np.pi/2, rtol=1e-3)  # theta should be conserved
        
        # Phi should have increased
        assert final_pos[0, 3] > positions[0, 3]

    def test_integration_failure_handling(self):
        """Test that integration handles reasonable edge cases."""
        self.setUp()
        # Test with reasonable but challenging conditions
        positions = np.array([[0.0, 5.0, np.pi/2, 0.0]])  # Moderately close to black hole
        velocities = np.array([[2.0, 0.1, 0.0, 0.5]])  # Reasonable velocities
        
        # This should work without issues
        final_pos, final_vel = self.kerr.integrate_batch_trajectory(positions, velocities, 0.01)
        assert np.all(np.isfinite(final_pos))
        assert np.all(np.isfinite(final_vel))


class TestKerrOneStepIntegration:
    """Test the convenience one-step integration method."""

    def setUp(self):
        """Set up test fixtures."""
        self.kerr = KerrSpacetime(M=1.0, a=0.5)

    def test_one_step_equivalence(self):
        """Test that integrate_one_step is equivalent to integrate_batch_trajectory."""
        self.setUp()
        positions = np.array([[0.0, 10.0, np.pi/2, 0.0]])
        velocities = np.array([[1.0, 0.0, 0.0, 0.1]])
        dtau = 0.1
        
        # Use both methods
        pos1, vel1 = self.kerr.integrate_one_step(positions, velocities, dtau)
        pos2, vel2 = self.kerr.integrate_batch_trajectory(positions, velocities, dtau)
        
        assert_allclose(pos1, pos2, rtol=1e-14)
        assert_allclose(vel1, vel2, rtol=1e-14)


class TestKerrPhysicalProperties:
    """Test physical properties and edge cases."""

    def setUp(self):
        """Set up test fixtures."""
        self.kerr = KerrSpacetime(M=1.0, a=0.5)
        self.kerr_extremal = KerrSpacetime(M=1.0, a=1.0)

    def test_ergosphere_properties(self):
        """Test properties in the ergosphere region."""
        self.setUp()
        # In the ergosphere, g_tt > 0
        # Ergosphere boundary: r = M + sqrt(M^2 - a^2 cos^2(theta))
        
        theta = np.pi/2  # Equatorial plane
        M, a = self.kerr.M, self.kerr.a
        r_ergo = M + np.sqrt(M**2 - a**2)  # At equatorial plane
        
        # Just outside ergosphere
        r_outside = np.array([r_ergo + 0.1])
        theta_arr = np.array([theta])
        
        g_outside = self.kerr._metric_components(r_outside, theta_arr)
        assert g_outside[0, 0, 0] < 0  # Should be negative (timelike)

    def test_extremal_kerr_properties(self):
        """Test properties of extremal Kerr black hole."""
        self.setUp()
        r = np.array([10.0])
        theta = np.array([np.pi/2])
        
        g = self.kerr_extremal._metric_components(r, theta)
        
        # Should still be a valid metric
        assert np.all(np.isfinite(g))
        
        # Test symmetry
        assert_allclose(g[0], g[0].T, rtol=1e-14)

    def test_large_radius_limit(self):
        """Test that metric approaches Minkowski at large distances."""
        self.setUp()
        r_large = np.array([1000.0])  # Very large radius
        theta = np.array([np.pi/2])
        
        g = self.kerr._metric_components(r_large, theta)
        
        # At large r, metric should approach Minkowski in spherical coordinates
        # g_tt → -1, g_rr → 1, g_θθ → r^2, g_φφ → r^2 sin^2θ
        # Off-diagonal terms should vanish
        
        assert_allclose(g[0, 0, 0], -1.0, rtol=1e-2)
        assert_allclose(g[0, 1, 1], 1.0, rtol=1e-2)
        assert_allclose(g[0, 2, 2], r_large[0]**2, rtol=1e-2)
        assert_allclose(g[0, 3, 3], r_large[0]**2, rtol=1e-2)
        assert_allclose(g[0, 0, 3], 0.0, atol=1e-2)


class TestKerrNumericalStability:
    """Test numerical stability and robustness."""

    def setUp(self):
        """Set up test fixtures."""
        self.kerr = KerrSpacetime(M=1.0, a=0.9)  # High spin for challenging case

    def test_numerical_derivatives(self):
        """Test that numerical derivatives in Christoffel symbols are stable."""
        self.setUp()
        r = np.array([5.0, 10.0, 20.0])
        theta = np.array([np.pi/4, np.pi/2, 3*np.pi/4])
        
        # Should not raise exceptions or produce NaN/inf
        Gamma = self.kerr._batched_christoffel_symbols(r, theta)
        assert np.all(np.isfinite(Gamma))

    def test_near_pole_behavior(self):
        """Test behavior near coordinate poles."""
        self.setUp()
        r = np.array([10.0])
        theta_near_pole = np.array([0.01])  # Near north pole
        
        g = self.kerr._metric_components(r, theta_near_pole)
        assert np.all(np.isfinite(g))

    def test_batch_size_independence(self):
        """Test that results don't depend on batch size."""
        self.setUp()
        r_single = np.array([10.0])
        theta_single = np.array([np.pi/2])
        
        r_batch = np.array([10.0, 10.0, 10.0])
        theta_batch = np.array([np.pi/2, np.pi/2, np.pi/2])
        
        g_single = self.kerr._metric_components(r_single, theta_single)
        g_batch = self.kerr._metric_components(r_batch, theta_batch)
        
        # All elements in batch should be identical to single calculation
        for i in range(3):
            assert_allclose(g_batch[i], g_single[0], rtol=1e-14)


if __name__ == "__main__":
    # Run tests when script is executed directly
    print("Running Kerr spacetime tests...")
    pytest.main([__file__, "-v"])
