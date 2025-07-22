#%%
"""
This notebook provides a hands-on explanation of Section 4.5 from the review
"The causal set approach to quantum gravity" by Sumati Surya. We will explore
how to calculate spacetime curvature (gravity) directly from the structure of
a causal set, leading to the famous Benincasa-Dowker action.
"""

#%%
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# Set a style for the plots
plt.style.use('seaborn-v0_8-whitegrid')

#%% md
# ### Introduction: Finding Gravity in the Dots
#
# In Einstein's General Relativity, gravity is not a force; it is the **curvature of spacetime**. The key mathematical object that describes this curvature at a point is the **Ricci Scalar, `R`**. The master equation for gravity, the Einstein-Hilbert action, is essentially the sum of `R` over all of spacetime.
#
# **The Challenge:** How can we possibly find a quantity like `R`, which describes smooth curvature, in a messy, random jumble of discrete dots?
#
# **The Solution:** We use a clever trick involving a discrete wave operator, `B`. By applying this operator to the most boring field imaginable (a field of all ones), the operator doesn't return zero as you'd expect. Instead, the non-zero result it spits out is directly proportional to the spacetime curvature `R`! This is the central insight of this section.

#%% md
# ### Step 1: Building the Tools (Recap from 4.4)
#
# To calculate curvature, we first need the tools to talk about "neighbors" in a causal set. The key idea is to classify every point in the past of a given element `e` based on how many other elements are "in between".
#
# #### Mathematical Foundation
#
# For a given element `e`, we define the set $L_k(e)$ as the set of all elements $e'$ in the past of $e$ such that the **causal interval** between them contains exactly `k` other elements.
#
# $$ L_k(e) = \{ e' \in C \mid e' \prec e \text{ and } |I(e', e)| = k \} $$
#
# *   $L_0(e)$ is the set of **links** to `e` (nearest neighbors).
# *   $L_1(e)$ is the set of **next-to-nearest neighbors**.
# *   ... and so on.
#
# We will need to count the number of elements in each of these sets. Let $N_k(e) = |L_k(e)|$.

#%%
# --- First, we need a causal set. We'll reuse our functions from before. ---

def get_causal_matrix(coords):
    """Computes the causal matrix for a set of 2D coordinates."""
    N = coords.shape[0]
    C = np.zeros((N, N), dtype=int)
    t_coords = coords[:, 0].reshape(N, 1)
    x_coords = coords[:, 1].reshape(N, 1)
    dt = t_coords - t_coords.T
    dx = x_coords - x_coords.T
    ds2 = dx**2 - dt**2
    causal_past_mask = (dt > 0) & (ds2 < 0)
    C[causal_past_mask] = 1
    return C

# --- Create a Causal Set ---
N_points = 100
coords = (np.random.rand(N_points, 2) - 0.5) * 2
diamond_coords = coords[np.sum(np.abs(coords), axis=1) < 1]
C = get_causal_matrix(diamond_coords)
N = diamond_coords.shape[0]

print(f"Created a causal set with {N} elements.")

#%%
def get_interval_sizes(causal_matrix):
    """
    For a causal matrix C, computes a new matrix where the entry (i, j)
    is the number of elements in the interval between i and j.
    """
    # The number of elements in the interval I(i, j) is the number of
    # elements k such that j < k < i.
    # This is equivalent to the number of paths of length 2 from j to i.
    # C[i, k] = 1 means k < i
    # C[k, j] = 1 means j < k
    # So, (C @ C)[i, j] = sum_k(C[i, k] * C[k, j]) is exactly what we need.
    return np.dot(causal_matrix, causal_matrix)

def get_Nk_counts_for_element(e_idx, causal_matrix, interval_sizes):
    """
    Calculates N_k(e) for k=0, 1, 2, 3 for a single element e.

    Args:
        e_idx (int): The index of the element 'e'.
        causal_matrix (np.ndarray): The full causal matrix.
        interval_sizes (np.ndarray): The matrix of interval sizes.

    Returns:
        dict: A dictionary with counts for N0, N1, N2, N3.
    """
    # Find all elements e' in the past of e
    past_indices = np.where(causal_matrix[e_idx, :] == 1)[0]

    # Get the sizes of the intervals |I(e, e')| for all e' in the past
    sizes_for_e = interval_sizes[e_idx, past_indices]

    # Count how many of these intervals have size 0, 1, 2, or 3
    counts = {
        0: np.sum(sizes_for_e == 0),
        1: np.sum(sizes_for_e == 1),
        2: np.sum(sizes_for_e == 2),
        3: np.sum(sizes_for_e == 3)
    }

    return counts

# Pre-calculate the interval sizes for the whole causet
interval_sizes_matrix = get_interval_sizes(C)

# Example: Get the neighbor counts for element 0
e_index_example = N - 1 # Pick the element with the largest time coord
Nk_example = get_Nk_counts_for_element(e_index_example, C, interval_sizes_matrix)

print(f"For element {e_index_example} (one of the latest points):")
print(f"  Number of nearest neighbors (links), N_0: {Nk_example[0]}")
print(f"  Number of next-nearest neighbors, N_1: {Nk_example[1]}")
print(f"  Number of N_2 neighbors: {Nk_example[2]}")
print(f"  Number of N_3 neighbors: {Nk_example[3]}")

#%% md
# ### Step 2: Calculating the Ricci Scalar `R`
#
# The paper presents the formula for the dimensionless discrete Ricci scalar `R` at an element `e` in 4-dimensions. It is a specific weighted sum of the number of its neighbors.
#
# #### Mathematical Foundation
#
# The formula (Eq. 33) is:
#
# $$ R(e) = \frac{4}{\sqrt{6}} \left[ 1 - N_0(e) + 9N_1(e) - 16N_2(e) + 8N_3(e) \right] $$
#
# *   **`R(e)`**: The Ricci scalar curvature at the specific dot `e`.
# *   **`N_k(e)`**: The number of neighbors of `e` at "distance" `k`, which we just calculated.
# *   **The Coefficients `(1, -1, 9, -16, 8)`**: These are not arbitrary! They are "magic numbers" derived from a detailed calculation in the continuum limit. They are precisely the weights needed for the non-local contributions to cancel out, leaving behind only the local curvature. The alternating signs are a hallmark of this cancellation.
# *   **$\frac{4}{\sqrt{6}}$**: This is a normalization constant for 4-dimensions.

#%%
def calculate_ricci_scalar_at_element(Nk_counts):
    """
    Calculates the discrete Ricci scalar for a single element e using its
    neighbor counts. (Formula for d=4)

    Args:
        Nk_counts (dict): A dictionary with keys 0,1,2,3 for N_k(e).

    Returns:
        float: The value of R(e).
    """
    # Magic coefficients from Eq. 33
    term = (1
            - Nk_counts[0]
            + 9 * Nk_counts[1]
            - 16 * Nk_counts[2]
            + 8 * Nk_counts[3])

    prefactor = 4 / np.sqrt(6)

    return prefactor * term

# Calculate R for our example element
R_example = calculate_ricci_scalar_at_element(Nk_example)
print(f"The discrete Ricci scalar curvature at element {e_index_example} is: R(e) = {R_example:.4f}")

#%% md
# ### Step 3: The Benincasa-Dowker (BD) Action
#
# Now that we can calculate the curvature `R` at any single point, we can find the total action for the entire universe. In General Relativity, this is the Einstein-Hilbert action. In Causal Set Theory, it's the **Benincasa-Dowker (BD) Action**.
#
# #### Mathematical Foundation
#
# The BD Action $S^{(4)}(C)$ for a 4D causal set `C` is simply the sum of the Ricci scalar over all elements in the set.
#
# $$ S^{(4)}(C) = \sum_{e \in C} R(e) $$
#
# This incredibly simple and elegant formula is the discrete, background-independent action for pure gravity. In the continuum-inspired dynamics of Chapter 6, this is the quantity that goes into the path integral ($e^{iS(c)}$) to determine the "quantum weight" of a given causal set universe.

#%%
def calculate_BD_action(causal_matrix):
    """
    Calculates the total Benincasa-Dowker action for a causal set.

    Args:
        causal_matrix (np.ndarray): The causal matrix for the set.

    Returns:
        float: The total action S(C).
    """
    total_action = 0
    num_elements = causal_matrix.shape[0]

    # Pre-calculate interval sizes to be efficient
    interval_sizes = get_interval_sizes(causal_matrix)

    # Loop over every element in the causal set
    for e_idx in range(num_elements):
        # Get neighbor counts for this element
        Nk_counts = get_Nk_counts_for_element(e_idx, causal_matrix, interval_sizes)

        # Calculate R(e) for this element
        R_e = calculate_ricci_scalar_at_element(Nk_counts)

        # Add it to the total sum
        total_action += R_e

    return total_action

# Calculate the total action for our toy universe
total_S_BD = calculate_BD_action(C)

print(f"The total Benincasa-Dowker Action for our {N}-element causal set is: S(C) = {total_S_BD:.4f}")

#%% md
# ### The "Smeared" Action (A More Robust Version)
#
# The paper notes that the formula for `R(e)` is very sensitive to Poisson fluctuations (the exact random placement of the dots). A more stable version uses a "smearing" function that averages over a range of interval sizes, controlled by a non-locality parameter `ε`.
#
# #### Mathematical Foundation
#
# Instead of a simple weighted sum of the $N_k$, the action is built from a more complicated operator $B_\epsilon$. The action takes the form:
#
# $$ S_\epsilon(C, \epsilon) = \sum_{e \in C} \left( -\phi(e) + \epsilon \sum_{e' \prec e} f(n(e', e), \epsilon) \phi(e') \right) $$
#
# *   **`ε`**: A new fundamental parameter of the theory, the "non-locality scale."
# *   **`n(e', e)`**: The number of elements in the interval $I(e', e)$.
# *   **`f(n, ε)`**: A specific "smearing function" (shown in Fig. 15) that replaces the simple integer coefficients. It smoothly weights the contributions from past points based on the interval size `n`.
#
# This version is much more stable and is used in the MCMC simulations discussed in Chapter 6. While we won't implement the full smearing function here, it's important to know that this more robust version exists.

#%% md
# ### Conclusion
#
# This notebook demonstrates the remarkable achievement of Section 4.5:
#
# 1.  We defined a way to count **"neighbors"** ($N_k$) in a causal set based purely on the causal order.
# 2.  Using a specific weighted sum of these neighbor counts, we defined a quantity `R(e)` at each dot that behaves exactly like the **Ricci scalar curvature** from General Relativity.
# 3.  Summing `R(e)` over all the dots gives us the **Benincasa-Dowker Action `S(C)`**, a complete, fundamental action for pure gravity on a discrete spacetime.
#
# This provides the final, crucial piece of the puzzle. We have successfully reconstructed the mathematical heart of Einstein's theory of gravity from the simple, discrete foundation of "Order + Number."