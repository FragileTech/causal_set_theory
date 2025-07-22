# %% [markdown]
"""
Of course! This is a fantastic way to learn. Section 5.1 is perfect for a hands-on demonstration because we can build a "toy universe" (a causal set) and then implement the Green functions directly.

Here is a Python notebook that explains the concepts from Section 5.1, using LaTeX for the math and providing functional Python implementations.
"""

# %%
# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# Set a style for the plots
plt.style.use('seaborn-v0_8-whitegrid')

# %% [markdown]
"""
### Introduction: What is a Green Function?

In physics, a **Green Function** is a tool that answers a simple question:

> If I poke spacetime at a point *A*, how much of a "ripple" do I feel at a point *B*?

It's a map of influence. It tells us how a disturbance (like the creation of a particle) propagates from a source point in the past to a destination point in the future. In standard physics, this is described by a smooth, continuous function.

**The Challenge:** How do we define this on a causal set, which is just a messy, random collection of discrete dots?

**The Solution:** We use the most fundamental information we have: the causal order itself. The influence can only travel from past dots to future dots. This notebook will show you how this simple idea is powerful enough to reconstruct the Green functions of continuum physics.
"""

# %% [markdown]
"""
### Step 1: Building Our Toy Universe (A 2D Causal Set)

Before we can define a Green function, we need a causal set to work with. We will create one by "sprinkling" `N` points randomly into a 2D Minkowski spacetime diamond and then determining their causal relationships.

#### Mathematical Foundation

A point $x_j = (t_j, x_j)$ is in the causal past of point $x_i = (t_i, x_i)$ if two conditions are met:
1.  **Time Ordering:** $t_i > t_j$ (The cause must precede the effect).
2.  **Timelike Separation:** The spacetime interval $\Delta s^2$ between them must be negative. The interval is defined as:
    $$ \Delta s^2 = (\Delta x)^2 - (\Delta t)^2 < 0 $$
    where $\Delta t = t_i - t_j$ and $\Delta x = x_i - x_j$. (We set the speed of light $c=1$).

We can store all of these relationships in a **Causal Matrix**, $C$.

$$
C_{ij} =
\begin{cases}
1 & \text{if event } j \text{ is in the causal past of event } i \\
0 & \text{otherwise}
\end{cases}
$$

This matrix **is** the complete "Order" information of our causal set.
"""

# %%
def get_causal_matrix(coords):
    """
    Computes the causal matrix for a set of 2D coordinates.

    Args:
        coords (np.ndarray): An (N, 2) array of (t, x) coordinates.

    Returns:
        np.ndarray: An (N, N) causal matrix C, where C[i, j] = 1
                    if j is in the causal past of i.
    """
    N = coords.shape[0]
    C = np.zeros((N, N), dtype=int)
    
    # Use broadcasting to efficiently calculate all differences
    # t_coords is (N, 1), t_coords.T is (1, N) -> dt is (N, N)
    t_coords = coords[:, 0].reshape(N, 1)
    x_coords = coords[:, 1].reshape(N, 1)
    
    dt = t_coords - t_coords.T
    dx = x_coords - x_coords.T
    
    # Calculate the squared interval for all pairs
    ds2 = dx**2 - dt**2
    
    # A point j is in the past of i if dt > 0 and the interval is timelike (ds2 < 0)
    causal_past_mask = (dt > 0) & (ds2 < 0)
    
    C[causal_past_mask] = 1
    
    return C

# %%
# --- Create and Visualize the Causal Set ---
N_points = 75
# Sprinkle points into a diamond shape |t| + |x| < 1
coords = (np.random.rand(N_points, 2) - 0.5) * 2
diamond_coords = coords[np.sum(np.abs(coords), axis=1) < 1]

# Get the causal matrix
C = get_causal_matrix(diamond_coords)
N = diamond_coords.shape[0]

# --- Visualization ---
fig, ax = plt.subplots(figsize=(10, 10))

# Plot the points
ax.scatter(diamond_coords[:, 0], diamond_coords[:, 1], c='crimson', zorder=2)

# Find the pairs of points that are causally related to draw lines
causally_related_indices = np.argwhere(C == 1)
lines = [[diamond_coords[j], diamond_coords[i]] for i, j in causally_related_indices]

line_collection = LineCollection(lines, colors='gray', linewidths=0.5, zorder=1, alpha=0.7)
ax.add_collection(line_collection)

ax.set_title(f'A 2D Causal Set with {N} Elements', fontsize=16)
ax.set_xlabel('Space (x)', fontsize=12)
ax.set_ylabel('Time (t)', fontsize=12)
ax.set_aspect('equal')
plt.show()

# %% [markdown]
"""
### Section 5.1a: The Green Function in 2 Dimensions

The paper explains that in 2D, the Green function is almost trivially related to the causal matrix we just built.

#### Mathematical Foundation

The causal set Green function, denoted $K_0^{(2)}$, is directly proportional to the causal matrix $C_0$ (which we called `C`).

$$ K_0^{(2)}(x, x') = \frac{1}{2} C_0(x, x') $$

This simple relationship comes from the fact that in 2D, the influence from a point is spread evenly across its entire future lightcone. The causal matrix captures this perfectly.
"""

# %%
def K0_2D(causal_matrix):
    """
    Calculates the 2D massless scalar Green function for a causal set.
    
    Args:
        causal_matrix (np.ndarray): The N x N causal matrix.
        
    Returns:
        np.ndarray: The N x N Green function matrix.
    """
    return 0.5 * causal_matrix

# Calculate the 2D Green function for our causal set
G_2D = K0_2D(C)

print("Shape of our Causal Matrix:", C.shape)
print("Shape of our 2D Green Function Matrix:", G_2D.shape)
print("\nExample: Influence of dot 0 on other dots (a row from the Green function matrix):")
# We transpose because G[i,j] means influence of j on i. 
# So influence *from* 0 is the 0-th column.
print(G_2D[:, 0])

# %% [markdown]
"""
### Section 5.1b: The Green Function in 4 Dimensions

In 4D (and other higher dimensions), the situation is more subtle. The influence
 of a disturbance is concentrated on the *edge* of the future lightcone, not spread 
 throughout its interior. The simple causal matrix is no longer a good enough tool.

We need a more refined object: the **Link Matrix**, $L_0$.

#### Mathematical Foundation

A causal relation $e' \prec e$ is called a **link** if there is no other
 element $z$ "in between" them, i.e., there is no $z$ such that $e' \prec z \prec e$. A link 
 represents an irreducible, direct step of causal influence.

We can find the link matrix from the causal matrix using matrix multiplication. 
The matrix $C^2 = C \cdot C$ counts the number of paths of length 2 between any two elements.
*   If $C_{ij} = 1$, there is a causal path from $j$ to $i$.
*   If $(C^2)_{ij} > 0$, there is at least one path of length 2 or more from $j$ to $i$.
*   Therefore, the relation from $j$ to $i$ is a **link** if and only if $C_{ij} = 1$ AND $(C^2)_{ij} = 0$.

The 4D Green function $K_0^{(4)}$ is then proportional to this link matrix:

$$ K_0^{(4)}(x, x') = \frac{1}{2\pi\sqrt{6}} L_0(x, x') $$

The constant prefactor $\frac{1}{2\pi\sqrt{6}}$ is a normalization constant that comes from a much more detailed calculation, ensuring the average over many sprinklings matches the continuum result. For us, it's just a constant.
"""

# %%
def get_link_matrix(causal_matrix):
    """
    Computes the link matrix from a causal matrix.

    Args:
        causal_matrix (np.ndarray): The N x N causal matrix.

    Returns:
        np.ndarray: The N x N link matrix L, where L[i, j] = 1
                    if j is a link to i.
    """
    # Matrix multiplication C @ C gives paths of length 2.
    C_squared = np.dot(causal_matrix, causal_matrix)
    
    # A link exists where a causal relation exists (C=1) but no
    # intermediate path exists (C_squared=0).
    is_a_link = (causal_matrix == 1) & (C_squared == 0)
    
    L = np.zeros_like(causal_matrix)
    L[is_a_link] = 1
    
    return L

def K0_4D(link_matrix):
    """
    Calculates the 4D massless scalar Green function for a causal set.
    
    Args:
        link_matrix (np.ndarray): The N x N link matrix.
        
    Returns:
        np.ndarray: The N x N Green function matrix.
    """
    prefactor = 1 / (2 * np.pi * np.sqrt(6))
    return prefactor * link_matrix

# Calculate the link matrix and the 4D Green function
L = get_link_matrix(C)
G_4D = K0_4D(L)

# %%
# --- Visualize the Links ---
fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(diamond_coords[:, 0], diamond_coords[:, 1], c='dodgerblue', zorder=2)
link_indices = np.argwhere(L == 1)
link_lines = [[diamond_coords[j], diamond_coords[i]] for i, j in link_indices]
link_line_collection = LineCollection(link_lines, colors='black', linewidths=0.8, zorder=1)
ax.add_collection(link_line_collection)
ax.set_title(f'Causal Set Links (Irreducible Relations)', fontsize=16)
ax.set_xlabel('Space (x)', fontsize=12)
ax.set_ylabel('Time (t)', fontsize=12)
ax.set_aspect('equal')
plt.show()

# %% [markdown]
"""
### Section 5.1c: Adding Mass - An Evocative Analogy

What happens if the particle is not massless? The paper describes a beautiful picture inspired by Feynman's path integral.

#### Conceptual Foundation

The Green function for a massive particle, $G_m$, can be thought of as a sum over all possible causal paths between two points. The massless Green function, $G_0$, only considers the most direct paths (links). Mass adds corrections.

The full massive propagator can be written as an infinite series:

$$ G_m = G_0 - m^2 (G_0 * G_0) + m^4 (G_0 * G_0 * G_0) - \dots $$

Here, `m` is the mass and `*` represents a "convolution" (sum over all intermediate points), which in our discrete case is just matrix multiplication.

The paper gives a lovely analogy for this:
*   A massive particle's path is a series of light-like **"hops"** (the links, represented by $G_0$).
*   At each point along the path, it can **"stop"** (the intermediate points), and at each stop, the mass `m` contributes to the calculation.
*   The full Green function is the sum over all possible paths with any number of stops.

This gives us a deep, microscopic picture of what mass *is* in this framework: it's an amplitude for a propagating particle to interact with the discrete points of spacetime itself.
"""

# %%
def G_massive_approx(G0, mass, order=2):
    """
    Calculates an approximation of the massive Green function.
    
    Args:
        G0 (np.ndarray): The massless Green function (e.g., G_2D or G_4D).
        mass (float): The mass of the particle.
        order (int): The number of terms to include in the series.
        
    Returns:
        np.ndarray: The approximate massive Green function matrix.
    """
    Gm = np.copy(G0)
    G0_power = np.copy(G0)
    
    for k in range(1, order + 1):
        # G0_power is G0*G0*...*G0 (k times)
        G0_power = np.dot(G0_power, G0)
        
        # Add the next term in the series
        term = ((-1)**k) * (mass**(2*k)) * G0_power
        Gm += term
        
    return Gm

# Let's calculate an approximation for a massive particle in our 2D universe
# We use G_2D as our massless propagator G0
mass_val = 0.1
G_massive = G_massive_approx(G_2D, mass=mass_val, order=5)

print("Shape of the approximate massive Green function:", G_massive.shape)
print("\nComparing massless vs massive influence FROM dot 0:")
print("Massless G_2D:", G_2D[:, 0])
print("Massive G_m: ", G_massive[:, 0].round(4))

# %% [markdown]
"""
### Conclusion

This notebook demonstrates the core ideas of Section 5.1:

1.  We can construct a "toy universe"—a causal set—from simple geometric principles.
2.  The fundamental **causal matrix `C`**, containing all the "Order" information, is directly proportional to the Green function in 2D.
3.  For higher dimensions, we need the more refined **link matrix `L`**, which captures direct, irreducible causal connections. This is also easily calculated from `C`.
4.  These tools provide a complete, microscopic description of how massless particles propagate on a discrete spacetime.
5.  The framework can be extended to massive particles through a beautiful analogy of "hops" and "stops," which corresponds to a matrix power series.

We have successfully rebuilt a key tool of quantum field theory from the ground up, using only the discrete dots and their causal relationships.
"""
