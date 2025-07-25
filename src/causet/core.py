"""
Core functions for Causal Set Theory

This module contains all the fundamental functions extracted from the Jupyter notebooks
for working with causal sets, including:
- Causal set generation (Poisson sprinkling)
- Green functions calculation
- Benincasa-Dowker action computation
- Sorkin-Johnston vacuum construction
- Entanglement entropy calculation
- Spectral dimension analysis

### Summary Table

The table below provides a direct link between the computational implementation in the Python code and the theoretical framework presented in the paper. This highlights how abstract physical and mathematical concepts are translated into concrete algorithms.

| Python Function/Class | Equation(s) in Paper | Description |
| :--- | :--- | :--- |
| **Causal Set Construction** | | |
| `get_causal_matrix` | (6), Definition on p. 11 | Implements the fundamental causal relation (`≺`). An element `j` is in the past of `i` if the spacetime interval is timelike or null (`ds² ≥ 0`) and the time coordinate is ordered (`Δt > 0`). This function builds the matrix representation of the poset from spacetime coordinates. |
| `get_link_matrix` | Definition on p. 16, used in (56), (59) | Computes the link matrix `L`, where `L(i, j) = 1` if `j` is a direct causal predecessor of `i` with no elements in between. This is equivalent to finding causal relations that are not bridged by any path of length two or more (`C²=0`). |
| `get_interval_sizes` | Definition on p. 11, used in (33), (34) | Calculates the number of elements in the causal interval `I[i, j]`, which is the set of elements causally between `j` and `i`. This is computationally equivalent to finding the number of paths of length 2 from `j` to `i`, given by `(C @ C)[i, j]`. |
| `sprinkle_minkowski_region` | (8), (9) | Implements the Poisson sprinkling process. The number of points `N` in a region of volume `V` is drawn from a Poisson distribution with mean `ρV`, where `ρ` is the density. |
| **Green Functions** | | |
| `K0_2D` | **(55)** | Calculates the exact dimensionless, massless retarded Green's function in 2D, which is directly proportional to the causal matrix `C₀`. |
| `K0_4D` | **(59)** | Calculates the approximation for the dimensionless, massless retarded Green's function in 4D, which is proportional to the link matrix `L₀`. |
| `G_massive_approx` | **(60), (63), (65)** | Implements the discrete version of the perturbative expansion for the massive Green's function `G_m` in terms of the massless one `G₀` and mass `m`. Matrix multiplication replaces the convolution operator (`*`). |
| **Benincasa-Dowker (BD) Action** | | |
| `calculate_ricci_scalar_at_element` | **(33)** | Computes the discrete Ricci scalar `R(e)` at a single element `e` using the number of its `k`-nearest neighbors (`N_k`) for `k=0,1,2,3`. |
| `calculate_BD_action` | **(34)** | Calculates the total Benincasa-Dowker action, which is the discrete analogue of the Einstein-Hilbert action. It is the sum of the Ricci scalars `R(e)` over all elements in the causal set. |
| **Quantum Field Theory & Sorkin-Johnston (SJ) Vacuum** | | |
| `get_iDelta_operator` | **(72)** | Constructs the Pauli-Jordan operator `iΔ`, which is defined as `i` times the difference between the retarded (`G_R`) and advanced (`G_A`) Green's functions. `G_R` is proportional to `C` and `G_A` is proportional to `Cᵀ`. |
| `get_SJ_Wightman_function` | **(75)** | Calculates the Sorkin-Johnston (SJ) vacuum Wightman function `W_SJ` by performing a spectral decomposition of the `iΔ` operator and summing over the contributions from its positive eigenvalues. |
| `calculate_ssee` | **(77)** | Calculates the Spacetime Entanglement Entropy (SSEE). This involves solving the generalized eigenvalue problem `W_A v = λ iΔ_A v` for a subregion `A` and then using the eigenvalues `λ` in the standard entropy formula for a fermionic system. |
| **Spectral Dimension** | | |
| `get_transition_matrix` | Concept in Sec. 5.4 | Creates a transition matrix `T` for a random walk on the undirected graph defined by the causal links (`L + Lᵀ`). This is a prerequisite for calculating the spectral dimension. |
| `calculate_spectral_dimension` | Concept in Sec. 5.4 | Estimates the spectral dimension `d_s` by analyzing the return probability of a random walk. It fits the scaling of the return probability `P(σ) ∼ σ^(-d_s/2)`, where `σ` is the number of steps. |
| **Curved Spacetime** | | |
| `DeSitterSpacetime.ds_to_minkowski` | Concept in Sec. 2 (p. 7-8) | Implements the conformal mapping from the de Sitter static patch to a region of Minkowski space, a key concept related to the Hawking-King-McCarthy-Malament (HKMM) theorem. |
| `CurvedSpacetime.is_causally_related` | Theorem 1 (HKMM) | The general principle that the causal structure determines the conformal geometry. For de Sitter space, this function uses the conformal map to check causality in the simpler embedding Minkowski space. |
| `CurvedSpacetime.sprinkle_and_build_causet` | (8), (9) | A generalized implementation of the Poisson sprinkling process for arbitrary curved spacetimes, using the region's volume and a method to generate points within it. |
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Protocol
import warnings

import numpy as np
from scipy.linalg import eig


# ============================================================================
# TYPE DEFINITIONS AND PROTOCOLS
# ============================================================================


class SpacetimeProtocol(Protocol):
    """Protocol defining the interface for spacetime objects."""

    dim: int

    def is_causally_related(self, p1: np.ndarray, p2: np.ndarray) -> bool:
        """Check if point p1 causally precedes point p2."""
        ...


# ============================================================================
# REGION DEFINITION DATACLASSES
# ============================================================================


@dataclass
class RegionDefinition(ABC):
    """
    Abstract base class for defining spacetime regions.

    This class provides the interface for all spacetime regions,
    ensuring consistent behavior across different region types.
    """

    @abstractmethod
    def get_volume(self, spacetime: SpacetimeProtocol) -> float:
        """
        Calculate the spacetime volume of this region in the given spacetime.

        Args:
            spacetime: The spacetime object containing geometry information

        Returns:
            The volume of the region
        """

    @abstractmethod
    def contains_point(self, point: np.ndarray, spacetime: SpacetimeProtocol) -> bool:
        """
        Check if a point is within this region in the given spacetime.

        Args:
            point: The point coordinates to check
            spacetime: The spacetime object

        Returns:
            True if the point is in the region, False otherwise
        """

    @abstractmethod
    def generate_bounding_box_point(self, spacetime: SpacetimeProtocol) -> np.ndarray:
        """
        Generate a random point in a bounding box that contains this region.

        Args:
            spacetime: The spacetime object

        Returns:
            Random point coordinates within the bounding box
        """

    def validate_spacetime_compatibility(self, spacetime: SpacetimeProtocol) -> bool:  # noqa: ARG002
        """
        Check if this region is compatible with the given spacetime.

        Args:
            spacetime: The spacetime object to check compatibility with

        Returns:
            True if compatible, False otherwise
        """
        return True  # Default implementation - can be overridden


# Minkowski spacetime regions


@dataclass
class MinkowskiDiamondRegion(RegionDefinition):
    """
    Causal diamond region in Minkowski space: |t| + |x| < radius

    This represents a causal diamond centered at the origin with the specified radius.
    """

    radius: float

    def __post_init__(self):
        """Validate parameters after initialization."""
        if self.radius <= 0:
            msg = "Radius must be positive"
            raise ValueError(msg)

    def get_volume(self, spacetime: SpacetimeProtocol) -> float:
        """Calculate the area of the diamond in 2D Minkowski space."""
        if spacetime.dim != 2:
            msg = "MinkowskiDiamondRegion only supports 2D spacetime"
            raise ValueError(msg)
        return 2 * self.radius**2

    def contains_point(self, point: np.ndarray, spacetime: SpacetimeProtocol) -> bool:  # noqa: ARG002
        """Check if point is inside the diamond."""
        if len(point) != 2:
            msg = "Point must be 2D for Minkowski diamond region"
            raise ValueError(msg)
        t, x = point
        return abs(t) + abs(x) < self.radius

    def generate_bounding_box_point(self, spacetime: SpacetimeProtocol) -> np.ndarray:  # noqa: ARG002
        """Generate a point uniformly in the square bounding box [-radius, radius]²."""
        return (np.random.rand(2) * 2 - 1) * self.radius

    def validate_spacetime_compatibility(self, spacetime: SpacetimeProtocol) -> bool:
        """Check compatibility with Minkowski spacetime."""
        return (
            spacetime.dim == 2
            and hasattr(spacetime, "__class__")
            and "Minkowski" in spacetime.__class__.__name__
        )


@dataclass
class MinkowskiRectangleRegion(RegionDefinition):
    """
    Rectangular region in Minkowski space: t_min < t < t_max, x_min < x < x_max

    This represents a rectangular spacetime region with specified bounds.
    """

    t_min: float
    t_max: float
    x_min: float
    x_max: float

    def __post_init__(self):
        """Validate parameters after initialization."""
        if self.t_max <= self.t_min:
            msg = "t_max must be greater than t_min"
            raise ValueError(msg)
        if self.x_max <= self.x_min:
            msg = "x_max must be greater than x_min"
            raise ValueError(msg)

    def get_volume(self, spacetime: SpacetimeProtocol) -> float:
        """Calculate the area of the rectangle in 2D Minkowski space."""
        if spacetime.dim != 2:
            msg = "MinkowskiRectangleRegion only supports 2D spacetime"
            raise ValueError(msg)
        return (self.t_max - self.t_min) * (self.x_max - self.x_min)

    def contains_point(self, point: np.ndarray, spacetime: SpacetimeProtocol) -> bool:  # noqa: ARG002
        """Check if point is inside the rectangle."""
        if len(point) != 2:
            msg = "Point must be 2D for Minkowski rectangle region"
            raise ValueError(msg)
        t, x = point
        return self.t_min < t < self.t_max and self.x_min < x < self.x_max

    def generate_bounding_box_point(self, spacetime: SpacetimeProtocol) -> np.ndarray:  # noqa: ARG002
        """Generate a point uniformly in the rectangle."""
        t = np.random.uniform(self.t_min, self.t_max)
        x = np.random.uniform(self.x_min, self.x_max)
        return np.array([t, x])

    def validate_spacetime_compatibility(self, spacetime: SpacetimeProtocol) -> bool:
        """Check compatibility with Minkowski spacetime."""
        return (
            spacetime.dim == 2
            and hasattr(spacetime, "__class__")
            and "Minkowski" in spacetime.__class__.__name__
        )


@dataclass
class MinkowskiLightconeRegion(RegionDefinition):
    """
    Future lightcone region in Minkowski space from point (t0, x0) up to t_max

    This represents the future lightcone emanating from a given point.
    """

    t0: float
    x0: float
    t_max: float

    def __post_init__(self):
        """Validate parameters after initialization."""
        if self.t_max <= self.t0:
            msg = "t_max must be greater than t0"
            raise ValueError(msg)

    def get_volume(self, spacetime: SpacetimeProtocol) -> float:
        """Calculate the area of the lightcone in 2D Minkowski space."""
        if spacetime.dim != 2:
            msg = "MinkowskiLightconeRegion only supports 2D spacetime"
            raise ValueError(msg)
        height = self.t_max - self.t0
        return height**2

    def contains_point(self, point: np.ndarray, spacetime: SpacetimeProtocol) -> bool:  # noqa: ARG002
        """Check if point is inside the lightcone."""
        if len(point) != 2:
            msg = "Point must be 2D for Minkowski lightcone region"
            raise ValueError(msg)
        t, x = point
        return self.t0 < t < self.t_max and abs(x - self.x0) < (t - self.t0)

    def generate_bounding_box_point(self, spacetime: SpacetimeProtocol) -> np.ndarray:  # noqa: ARG002
        """Generate a point uniformly in the triangular bounding box."""
        height = self.t_max - self.t0
        t = np.random.uniform(self.t0, self.t_max)
        x = np.random.uniform(self.x0 - height, self.x0 + height)
        return np.array([t, x])

    def validate_spacetime_compatibility(self, spacetime: SpacetimeProtocol) -> bool:
        """Check compatibility with Minkowski spacetime."""
        return (
            spacetime.dim == 2
            and hasattr(spacetime, "__class__")
            and "Minkowski" in spacetime.__class__.__name__
        )


# de Sitter spacetime regions


@dataclass
class DeSitterStaticPatchRegion(RegionDefinition):
    """
    Full static patch region in de Sitter space

    This represents the entire static patch of de Sitter spacetime,
    bounded by the horizons.
    """

    def get_volume(self, spacetime: SpacetimeProtocol) -> float:
        """
        Calculate the volume of the static patch.

        Volume element is dV = α*cosh(τ/α) dτ dχ
        Static patch boundaries are τ/α ± χ = ±π
        This integral evaluates to α² * π²
        """
        if spacetime.dim != 2:
            msg = "DeSitterStaticPatchRegion only supports 2D spacetime"
            raise ValueError(msg)
        if not hasattr(spacetime, "alpha"):
            msg = "Spacetime must have 'alpha' attribute for de Sitter calculations"
            raise ValueError(msg)
        return spacetime.alpha**2 * np.pi**2

    def contains_point(self, point: np.ndarray, spacetime: SpacetimeProtocol) -> bool:  # noqa: ARG002
        """
        Check if point is in the static patch.

        Our generation method always produces valid points in the static patch,
        so this always returns True.
        """
        if len(point) != 2:
            msg = "Point must be 2D for de Sitter static patch region"
            raise ValueError(msg)
        return True

    def generate_bounding_box_point(self, spacetime: SpacetimeProtocol) -> np.ndarray:
        """
        Generate a point uniformly in the static patch.

        Generate in (u,v) = (tau/α+chi, tau/α-chi) coordinates
        where -π < u < π and -π < v < π
        """
        if not hasattr(spacetime, "alpha"):
            msg = "Spacetime must have 'alpha' attribute"
            raise ValueError(msg)

        u = np.random.uniform(-np.pi, np.pi)
        v = np.random.uniform(-np.pi, np.pi)

        tau_over_alpha = (u + v) / 2
        chi = (u - v) / 2

        return np.array([tau_over_alpha * spacetime.alpha, chi])

    def validate_spacetime_compatibility(self, spacetime: SpacetimeProtocol) -> bool:
        """Check compatibility with de Sitter spacetime."""
        return (
            spacetime.dim == 2
            and hasattr(spacetime, "alpha")
            and hasattr(spacetime, "__class__")
            and "DeSitter" in spacetime.__class__.__name__
        )


@dataclass
class DeSitterDiamondRegion(RegionDefinition):
    """
    Diamond region in de Sitter space: |tau/alpha + chi| < radius, |tau/alpha - chi| < radius

    This represents a causal diamond in de Sitter spacetime with specified radius.
    """

    radius: float

    def __post_init__(self):
        """Validate parameters after initialization."""
        if self.radius <= 0:
            msg = "Radius must be positive"
            raise ValueError(msg)
        if self.radius > np.pi:
            warnings.warn("Radius larger than π may exceed static patch boundaries")

    def get_volume(self, spacetime: SpacetimeProtocol) -> float:
        """Calculate the area of the diamond region in de Sitter space."""
        if spacetime.dim != 2:
            msg = "DeSitterDiamondRegion only supports 2D spacetime"
            raise ValueError(msg)
        if not hasattr(spacetime, "alpha"):
            msg = "Spacetime must have 'alpha' attribute for de Sitter calculations"
            raise ValueError(msg)
        # Volume scales as radius² for the diamond region
        return spacetime.alpha**2 * self.radius**2

    def contains_point(self, point: np.ndarray, spacetime: SpacetimeProtocol) -> bool:
        """Check if point is inside the diamond region."""
        if len(point) != 2:
            msg = "Point must be 2D for de Sitter diamond region"
            raise ValueError(msg)
        if not hasattr(spacetime, "alpha"):
            msg = "Spacetime must have 'alpha' attribute"
            raise ValueError(msg)

        tau, chi = point
        tau_over_alpha = tau / spacetime.alpha
        return abs(tau_over_alpha + chi) < self.radius and abs(tau_over_alpha - chi) < self.radius

    def generate_bounding_box_point(self, spacetime: SpacetimeProtocol) -> np.ndarray:
        """Generate a point uniformly in the diamond region."""
        if not hasattr(spacetime, "alpha"):
            msg = "Spacetime must have 'alpha' attribute"
            raise ValueError(msg)

        u = np.random.uniform(-self.radius, self.radius)
        v = np.random.uniform(-self.radius, self.radius)

        tau_over_alpha = (u + v) / 2
        chi = (u - v) / 2

        return np.array([tau_over_alpha * spacetime.alpha, chi])

    def validate_spacetime_compatibility(self, spacetime: SpacetimeProtocol) -> bool:
        """Check compatibility with de Sitter spacetime."""
        return (
            spacetime.dim == 2
            and hasattr(spacetime, "alpha")
            and hasattr(spacetime, "__class__")
            and "DeSitter" in spacetime.__class__.__name__
        )


# Type alias for any region definition
RegionType = (
    MinkowskiDiamondRegion
    | MinkowskiRectangleRegion
    | MinkowskiLightconeRegion
    | DeSitterStaticPatchRegion
    | DeSitterDiamondRegion
)
# ============================================================================
# BASIC CAUSAL SET CONSTRUCTION
# ============================================================================


def sprinkle_minkowski_diamond(radius: float, density: float = 1.0) -> np.ndarray:
    """
    Performs a Poisson sprinkling into a 2D Minkowski diamond |t|+|x| < radius.

    Args:
        radius (float): The radius of the causal diamond.
        density (float): The sprinkling density.

    Returns:
        np.ndarray: An (N, 2) array of sprinkled (t, x) coordinates.

    Raises:
        ValueError: If radius is non-positive or density is non-positive.
    """
    if radius <= 0:
        msg = "Radius must be positive"
        raise ValueError(msg)
    if density <= 0:
        msg = "Density must be positive"
        raise ValueError(msg)

    # Use the new MinkowskiDiamondRegion dataclass for consistency
    region = MinkowskiDiamondRegion(radius=radius)
    return sprinkle_minkowski_region(region, density)


def sprinkle_minkowski_region(region: RegionDefinition, density: float = 1.0) -> np.ndarray:
    """
    Performs a Poisson sprinkling into a specified 2D Minkowski region.

    Args:
        region (RegionDefinition): The region definition object.
        density (float): The sprinkling density.

    Returns:
        np.ndarray: An (N, 2) array of sprinkled (t, x) coordinates.

    Raises:
        ValueError: If density is non-positive or region is incompatible.

    Examples:
        # Diamond region
        diamond = MinkowskiDiamondRegion(radius=2.0)
        coords = sprinkle_minkowski_region(diamond, density=3.0)

        # Rectangular region
        rect = MinkowskiRectangleRegion(t_min=0, t_max=5, x_min=-2, x_max=2)
        coords = sprinkle_minkowski_region(rect, density=2.0)

        # Future lightcone
        lightcone = MinkowskiLightconeRegion(t0=0, x0=0, t_max=3)
        coords = sprinkle_minkowski_region(lightcone, density=4.0)
        
    References:
        Eq. (8) and (9): Implements the Poisson sprinkling process where 
        N ~ Poisson(ρV) for density ρ and region volume V.
    """
    if density <= 0:
        msg = "Density must be positive"
        raise ValueError(msg)

    spacetime = MinkowskiSpacetime()

    # Validate region compatibility
    if not region.validate_spacetime_compatibility(spacetime):
        msg = f"Region {type(region).__name__} is not compatible with MinkowskiSpacetime"
        raise ValueError(msg)

    coords, _ = spacetime.sprinkle_and_build_causet(region, density)
    return coords


def get_causal_matrix(coords: np.ndarray) -> np.ndarray:
    """
    Computes the causal matrix for a set of 2D coordinates.

    Args:
        coords (np.ndarray): An (N, 2) array of (t, x) coordinates.

    Returns:
        np.ndarray: An (N, N) causal matrix C, where C[i, j] = 1
                    if j is in the causal past of i.

    Raises:
        ValueError: If coords has wrong shape or contains invalid values.
        
    References:
        Eq. (6) and Definition on p. 11: Implements the fundamental causal 
        relation (≺) for building the matrix representation of the poset.
    """
    if coords.size == 0:
        return np.array([]).reshape(0, 0)

    if coords.ndim != 2 or coords.shape[1] != 2:
        msg = "Coordinates must be an (N, 2) array"
        raise ValueError(msg)

    if not np.isfinite(coords).all():
        msg = "Coordinates must contain only finite values"
        raise ValueError(msg)

    N = coords.shape[0]

    # Use broadcasting to efficiently calculate all differences.
    # coords[:,0] is a 1D array of t-values. Reshape to (N,1) for broadcasting.
    t_coords = coords[:, 0].reshape(N, 1)
    x_coords = coords[:, 1].reshape(N, 1)

    # dt[i,j] = t_i - t_j
    dt = t_coords - t_coords.T
    dx = x_coords - x_coords.T

    # Calculate the squared interval for all pairs
    interval_sq = dt**2 - dx**2

    # A point j is in the past of i if dt > 0 and the interval is timelike/null.
    causal_past_mask = (dt > 0) & (interval_sq >= 0)

    causal_matrix = np.zeros((N, N), dtype=int)
    causal_matrix[causal_past_mask] = 1

    return causal_matrix


def get_link_matrix(causal_matrix: np.ndarray) -> np.ndarray:
    """
    Computes the link matrix from a causal matrix.
    
    The link matrix L is defined as:
    L(x,x') = 1 if x' ≺ x is a link, 0 otherwise
    
    A link between x' and x means that x' ≺ x but there is no element z
    such that x' ≺ z ≺ x. This implements the discrete analogue of the
    lightcone structure used in Green function calculations.

    Args:
        causal_matrix (np.ndarray): The N x N causal matrix.

    Returns:
        np.ndarray: The N x N link matrix L, where L[i, j] = 1
                    if j is a link to i.

    Raises:
        ValueError: If causal_matrix has wrong shape or invalid values.
        
    References:
        Definition on p. 16, used in Eq. (56) and (59): Computes the link 
        matrix L for discrete analogue of the lightcone structure used in 
        Green function calculations.
    """
    if causal_matrix.size == 0:
        return np.array([]).reshape(0, 0)

    if causal_matrix.ndim != 2 or causal_matrix.shape[0] != causal_matrix.shape[1]:
        msg = "Causal matrix must be square"
        raise ValueError(msg)

    if not np.all((causal_matrix == 0) | (causal_matrix == 1)):
        msg = "Causal matrix must contain only 0s and 1s"
        raise ValueError(msg)

    # Matrix multiplication C @ C gives paths of length 2.
    C_squared = np.dot(causal_matrix, causal_matrix)

    # A link exists where a causal relation exists (C=1) but no
    # intermediate path exists (C_squared=0).
    is_a_link = (causal_matrix == 1) & (C_squared == 0)

    L = np.zeros_like(causal_matrix)
    L[is_a_link] = 1

    return L


def get_interval_sizes(causal_matrix: np.ndarray) -> np.ndarray:
    """
    For a causal matrix C, computes a new matrix where the entry (i, j)
    is the number of elements in the interval between i and j.

    Args:
        causal_matrix (np.ndarray): The N x N causal matrix.

    Returns:
        np.ndarray: The N x N interval size matrix.

    Raises:
        ValueError: If causal_matrix has wrong shape or invalid values.
        
    References:
        Definition on p. 11, used in Eq. (33) and (34): Calculates the number 
        of elements in causal intervals I[i,j] for Benincasa-Dowker action.
    """
    if causal_matrix.size == 0:
        return np.array([]).reshape(0, 0)

    if causal_matrix.ndim != 2 or causal_matrix.shape[0] != causal_matrix.shape[1]:
        msg = "Causal matrix must be square"
        raise ValueError(msg)

    if not np.all((causal_matrix == 0) | (causal_matrix == 1)):
        msg = "Causal matrix must contain only 0s and 1s"
        raise ValueError(msg)

    # The number of elements in the interval I(i, j) is the number of
    # elements k such that j < k < i.
    # This is equivalent to the number of paths of length 2 from j to i.
    # C[i, k] = 1 means k < i
    # C[k, j] = 1 means j < k
    # So, (C @ C)[i, j] = sum_k(C[i, k] * C[k, j]) is exactly what we need.
    return np.dot(causal_matrix, causal_matrix)


# ============================================================================
# GREEN FUNCTIONS
# ============================================================================


def K0_2D(causal_matrix: np.ndarray) -> np.ndarray:
    """
    Calculates the 2D massless scalar Green function for a causal set.
    
    Implements the discrete retarded Green function in 2D:
    K₀⁽²⁾(x,x') = (1/2) C₀(x,x')
    
    where C₀ is the causal matrix. This gives the exact massless retarded 
    Green function for causal sets sprinkled in 2D Minkowski spacetime.

    Args:
        causal_matrix (np.ndarray): The N x N causal matrix C₀.

    Returns:
        np.ndarray: The N x N Green function matrix K₀⁽²⁾.
        
    References:
        Eq. (55): K₀⁽²⁾(x,x') = (1/2)C₀(x,x') - Exact dimensionless, massless 
        retarded Green's function in 2D.
    """
    return 0.5 * causal_matrix


def K0_4D(link_matrix: np.ndarray) -> np.ndarray:
    """
    Calculates the 4D massless scalar Green function for a causal set.
    
    Implements the discrete retarded Green function in 4D:
    K₀⁽⁴⁾(x,x') = (1/(2π√6)) L₀(x,x')
    
    where L₀ is the link matrix. In the continuum limit ρ_c → ∞, this
    approximates the 4D massless retarded Green function.

    Args:
        link_matrix (np.ndarray): The N x N link matrix L₀.

    Returns:
        np.ndarray: The N x N Green function matrix K₀⁽⁴⁾.
        
    References:
        Eq. (59): K₀⁽⁴⁾(x,x') = (1/(2π√6)) L₀(x,x') - Approximation for the 
        dimensionless, massless retarded Green's function in 4D.
    """
    prefactor = 1 / (2 * np.pi * np.sqrt(6))
    return prefactor * link_matrix


def G_massive_approx(G0: np.ndarray, mass: float, order: int = 2) -> np.ndarray:
    """
    Calculates an approximation of the massive Green function.

    Args:
        G0 (np.ndarray): The massless Green function (e.g., from K0_2D or K0_4D).
        mass (float): The mass of the particle.
        order (int): The number of terms to include in the series.

    Returns:
        np.ndarray: The approximate massive Green function matrix.
        
    References:
        Eq. (60), (63), and (65): Discrete version of the perturbative expansion 
        for the massive Green's function G_m in terms of massless G₀ and mass m.
    """
    if G0.shape[0] == 0:
        return np.array([]).reshape(0, 0)

    Gm = np.copy(G0)
    G0_power = np.copy(G0)

    for k in range(1, order + 1):
        # G0_power is G0*G0*...*G0 (k times)
        G0_power = np.dot(G0_power, G0)

        # Add the next term in the series
        term = ((-1) ** k) * (mass ** (2 * k)) * G0_power
        Gm += term

    return Gm


# ============================================================================
# BENINCASA-DOWKER ACTION
# ============================================================================


def get_Nk_counts_for_element(
    e_idx: int, causal_matrix: np.ndarray, interval_sizes: np.ndarray
) -> dict:
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

    if len(past_indices) == 0:
        return {0: 0, 1: 0, 2: 0, 3: 0}

    # Get the sizes of the intervals |I(e, e')| for all e' in the past
    sizes_for_e = interval_sizes[e_idx, past_indices]

    # Count how many of these intervals have size 0, 1, 2, or 3
    return {
        0: np.sum(sizes_for_e == 0),
        1: np.sum(sizes_for_e == 1),
        2: np.sum(sizes_for_e == 2),
        3: np.sum(sizes_for_e == 3),
    }


def calculate_ricci_scalar_at_element(Nk_counts: dict) -> float:
    """
    Calculates the discrete Ricci scalar for a single element e using its
    neighbor counts. (Formula for d=4)
    
    Implements the Benincasa-Dowker discrete Ricci scalar:
    R(e) = (4/√6) × [1 - N₀(e) + 9N₁(e) - 16N₂(e) + 8N₃(e)]
    
    where Nₖ(e) is the number of k-element intervals that contain element e.
    This is the discrete analogue of the continuum Ricci scalar curvature.

    Args:
        Nk_counts (dict): A dictionary with keys 0,1,2,3 for N_k(e).

    Returns:
        float: The value of R(e).
        
    References:
        Eq. (33): R(e) = (4/√6) × [1 - N₀(e) + 9N₁(e) - 16N₂(e) + 8N₃(e)]
        Computes the discrete Ricci scalar at a single element using its 
        k-nearest neighbors for the Benincasa-Dowker action formula.
    """
    # Magic coefficients from Eq. 33
    term = 1 - Nk_counts[0] + 9 * Nk_counts[1] - 16 * Nk_counts[2] + 8 * Nk_counts[3]

    prefactor = 4 / np.sqrt(6)

    return prefactor * term


def calculate_BD_action(causal_matrix: np.ndarray) -> float:
    """
    Calculates the total Benincasa-Dowker action for a causal set.
    
    Implements the discrete Einstein-Hilbert action:
    S(C) = Σₑ R(e)
    
    where R(e) is the discrete Ricci scalar at element e. This action
    provides a discrete analogue of the continuum Einstein-Hilbert action
    and is used in causal set quantum gravity path integrals.

    Args:
        causal_matrix (np.ndarray): The causal matrix for the set.

    Returns:
        float: The total action S(C).
        
    References:
        Eq. (34): S(C) = Σₑ R(e) - Total Benincasa-Dowker action as the 
        discrete analogue of the Einstein-Hilbert action for causal set 
        quantum gravity path integrals.
    """
    if causal_matrix.shape[0] == 0:
        return 0.0

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


# ============================================================================
# SORKIN-JOHNSTON VACUUM AND QUANTUM FIELD THEORY
# ============================================================================


def get_iDelta_operator(causal_matrix: np.ndarray) -> np.ndarray:
    """
    Constructs the Pauli-Jordan operator iΔ for a massless 2D causal set.
    
    Implements the discrete Pauli-Jordan function:
    iΔ = i(G_R - G_A)
    
    where G_R and G_A are the retarded and advanced Green functions.
    In 2D: G_R ∝ C and G_A ∝ C^T, where C is the causal matrix.
    
    This operator is fundamental for the Sorkin-Johnston vacuum construction
    and provides the commutation relations for quantum field theory on causal sets.

    Args:
        causal_matrix (np.ndarray): The causal matrix C.

    Returns:
        np.ndarray: The Pauli-Jordan operator iΔ.
        
    References:
        Eq. (72): iΔ = i(G_R - G_A) - Constructs the Pauli-Jordan operator 
        for the Sorkin-Johnston vacuum construction and quantum field theory 
        commutation relations on causal sets.
    """
    if causal_matrix.shape[0] == 0:
        return np.array([]).reshape(0, 0)

    # G_R is influence from past to future (j -> i), proportional to C[i,j]
    G_R = 0.5 * causal_matrix
    # G_A is influence from future to past (j -> i), proportional to C[j,i]
    G_A = 0.5 * causal_matrix.T

    Delta = G_R - G_A
    return 1j * Delta


def get_SJ_Wightman_function(iDelta_op: np.ndarray) -> np.ndarray:
    """
    Calculates the Sorkin-Johnston vacuum Wightman function from the iΔ operator.
    
    Implements the SJ vacuum construction:
    W_SJ(x,x') = Σ_{λₖ>0} λₖ |ψₖ⟩⟨ψₖ|
    
    where λₖ and |ψₖ⟩ are the positive eigenvalues and eigenvectors of iΔ.
    This gives the two-point function ⟨0|φ(x)φ(x')|0⟩ in the SJ vacuum state.
    
    The construction follows from requiring that the vacuum minimizes the 
    renormalized expectation value of ∇² while respecting the discrete 
    causal structure.

    Args:
        iDelta_op (np.ndarray): The iΔ (Pauli-Jordan) operator matrix.

    Returns:
        np.ndarray: The N × N Wightman function matrix W_SJ.
        
    References:
        Eq. (75): W_SJ(x,x') = Σ_{λₖ>0} λₖ |ψₖ⟩⟨ψₖ| - Sorkin-Johnston vacuum 
        Wightman function from spectral decomposition of the iΔ operator.
    """
    if iDelta_op.shape[0] == 0:
        return np.array([]).reshape(0, 0)

    # Step 1: Find eigenvalues and eigenvectors of the Hermitian operator iΔ
    eigenvalues, eigenvectors = np.linalg.eigh(iDelta_op)

    # Step 2: Construct the Wightman function by summing over positive eigenvalues
    W_SJ = np.zeros_like(iDelta_op, dtype=complex)
    for i in range(len(eigenvalues)):
        if eigenvalues[i] > 1e-9:  # Use a small tolerance for zero
            lambda_k = eigenvalues[i]
            v_k = eigenvectors[:, i]
            # Add the contribution from this positive mode: λ * v * v†
            W_SJ += lambda_k * np.outer(v_k, np.conj(v_k))

    return W_SJ


def calculate_ssee(W_SJ: np.ndarray, iDelta: np.ndarray, region_indices: list) -> float:
    """
    Calculates the Spacetime Entanglement Entropy (SSEE) for a subregion.
    
    Implements the generalized eigenvalue problem:
    W_A |ψₖ⟩ = λₖ iΔ_A |ψₖ⟩
    
    followed by the entropy calculation:
    S = -Σₖ [λₖ ln(λₖ) + (1-λₖ) ln(1-λₖ)]
    
    where W_A and iΔ_A are the restrictions of the Wightman function and 
    Pauli-Jordan operator to the subregion A. The eigenvalues λₖ represent
    the "occupation probabilities" of the modes in the subregion.
    
    This generalizes the usual entanglement entropy to quantum field theory
    on causal sets, providing a spacetime-based measure of entanglement.

    Args:
        W_SJ (np.ndarray): The full Wightman function for the causal set.
        iDelta (np.ndarray): The full iΔ operator for the causal set.
        region_indices (list or np.ndarray): A list of indices of the points
                                             that form the subregion 'A'.

    Returns:
        float: The entanglement entropy S.
        
    References:
        Eq. (77): S = -Σₖ [λₖ ln(λₖ) + (1-λₖ) ln(1-λₖ)] - Spacetime 
        Entanglement Entropy from generalized eigenvalue problem for 
        subregions in causal set quantum field theory.
    """
    if len(region_indices) == 0:
        return 0.0

    # Step 1: Restrict the operators to the subregion A
    W_A = W_SJ[np.ix_(region_indices, region_indices)]
    iDelta_A = iDelta[np.ix_(region_indices, region_indices)]

    # Step 2: Solve the generalized eigenvalue problem W_A * v = λ * (iΔ_A) * v
    try:
        # The 'eig' function can handle generalized eigenvalue problems
        gen_eigenvalues, _ = eig(W_A, b=iDelta_A)
    except np.linalg.LinAlgError:
        return np.nan

    # Step 3: Calculate the entropy from the eigenvalues
    entropy = 0.0
    for lam in gen_eigenvalues:
        # Physical eigenvalues for this problem should be real and between 0 and 1.
        if 0 < np.real(lam) < 1 and np.isclose(np.imag(lam), 0):
            lambda_k = np.real(lam)
            # Formula: S = -Sum[p*ln(p) + (1-p)*ln(1-p)]
            entropy -= lambda_k * np.log(lambda_k) + (1 - lambda_k) * np.log(1 - lambda_k)

    return entropy


# ============================================================================
# SPECTRAL DIMENSION
# ============================================================================


def get_transition_matrix(link_matrix: np.ndarray) -> np.ndarray:
    """
    Creates a transition matrix for a random walk on the causal set graph.
    
    Constructs the symmetric transition matrix:
    T_{ij} = A_{ij} / deg(i)
    
    where A = L + L^T is the undirected adjacency matrix from the link matrix,
    and deg(i) is the degree of vertex i. This gives equal probability of 
    transitioning to any neighboring vertex.
    
    The resulting matrix is used in spectral dimension calculations via
    the analysis of return probabilities in discrete random walks.

    Args:
        link_matrix (np.ndarray): The N × N link matrix L of the causal set.

    Returns:
        np.ndarray: The N × N transition matrix T for the random walk.
        
    References:
        Concept in Sec. 5.4: Creates transition matrix T for random walk 
        on undirected graph defined by causal links for spectral dimension 
        calculations.
    """
    if link_matrix.shape[0] == 0:
        return np.array([]).reshape(0, 0)

    # Create an undirected adjacency matrix A = L + L^T
    adjacency_matrix = link_matrix + link_matrix.T
    degree = np.sum(adjacency_matrix, axis=0)
    degree[degree == 0] = 1

    return adjacency_matrix / degree


def calculate_spectral_dimension(link_matrix: np.ndarray, max_steps: int = 50) -> float:
    """
    Calculates the spectral dimension via discrete random walk analysis.
    
    Implements the spectral dimension calculation:
    d_s = -2 × d(ln⟨P_σ⟩)/d(ln σ)
    
    where ⟨P_σ⟩ is the average return probability after σ steps of a random walk.
    The spectral dimension characterizes the effective dimensionality of the 
    discrete geometry as probed by diffusion processes.
    
    The calculation involves:
    1. Constructing the transition matrix T from the link structure
    2. Computing return probabilities ⟨P_σ⟩ = Tr(T^σ)/N for various σ
    3. Fitting the scaling ⟨P_σ⟩ ∼ σ^(-d_s/2) to extract d_s

    Args:
        link_matrix (np.ndarray): The N × N link matrix of the causal set.
        max_steps (int): Maximum number of random walk steps to simulate.

    Returns:
        float: The calculated spectral dimension d_s.
        
    References:
        Concept in Sec. 5.4: d_s = -2 × d(ln⟨P_σ⟩)/d(ln σ) - Estimates 
        spectral dimension by analyzing return probability scaling P(σ) ∼ σ^(-d_s/2) 
        in random walks on the discrete geometry.
    """
    if link_matrix.shape[0] == 0:
        return np.nan

    # Step 1: Get the transition matrix
    T = get_transition_matrix(link_matrix)

    # Step 2: Calculate the Heat Kernel
    avg_return_prob = []
    steps = np.arange(1, max_steps + 1)

    T_power_sigma = np.copy(T)
    for sigma in steps:
        if sigma > 1:
            T_power_sigma = np.dot(T_power_sigma, T)

        # The return probability is the diagonal of the matrix
        return_probabilities = np.diag(T_power_sigma)
        avg_return_prob.append(np.mean(return_probabilities))

    avg_return_prob = np.array(avg_return_prob)

    # Step 3: Fit a line to the log-log plot to find the slope
    # Fit for steps in the middle range to get stable results
    fit_start = max(1, max_steps // 4)
    fit_end = min(len(steps), 3 * max_steps // 4)

    if fit_end <= fit_start:
        return np.nan

    log_steps = np.log(steps[fit_start:fit_end])
    log_return_prob = np.log(avg_return_prob[fit_start:fit_end])

    # Filter out infinite or NaN values
    valid_mask = np.isfinite(log_steps) & np.isfinite(log_return_prob)
    if np.sum(valid_mask) < 2:
        return np.nan

    log_steps = log_steps[valid_mask]
    log_return_prob = log_return_prob[valid_mask]

    # Fit a linear polynomial (degree 1)
    slope, _ = np.polyfit(log_steps, log_return_prob, 1)

    # The spectral dimension d_s = -2 * slope
    return -2 * slope


# ============================================================================
# CURVED SPACETIME FRAMEWORK
# ============================================================================


class CurvedSpacetime(ABC):
    """
    Abstract base class for curved spacetime geometries in causal set theory.
    
    This class provides the framework for implementing causal set constructions
    in arbitrary curved spacetimes. Concrete implementations must specify:
    
    1. The causal structure via is_causally_related()
    2. The conformal factor and metric properties
    3. Region-specific sprinkling algorithms
    
    The class supports the general causal set sprinkling procedure:
    - Poisson process with density ρ in the spacetime volume
    - Causal matrix construction from the spacetime geometry  
    - Integration with RegionDefinition dataclasses
    
    Subclasses implement specific geometries like Minkowski, de Sitter,
    anti-de Sitter, or more general curved backgrounds.

    Attributes:
        dim (int): The spacetime dimension.
        
    References:
        General framework for causal sets in curved spacetime following
        the covariant approach of Sorkin and others. The sprinkling density
        is determined by the spacetime volume element √|g| d^n x.
    """

    def __init__(self, dim: int):
        self.dim = dim

    @abstractmethod
    def is_causally_related(self, p1: np.ndarray, p2: np.ndarray) -> bool:
        """
        Checks for a causal relationship between two spacetime points.
        
        Determines whether point p1 causally precedes point p2 according to
        the spacetime metric. This implements the fundamental causal structure
        that defines the partial ordering in causal set theory.
        
        For a general metric g_μν, the causal relation requires:
        1. The connecting vector to be timelike or null: g_μν Δx^μ Δx^ν ≥ 0
        2. The time ordering to be correct: Δt > 0
        
        Args:
            p1, p2 (np.ndarray): Spacetime coordinate arrays.
            
        Returns:
            bool: True if p1 ≺ p2 (p1 causally precedes p2).
        """

    def sprinkle_and_build_causet(self, region: RegionDefinition, density: float = 1.0):
        """
        Performs a Poisson sprinkling and constructs the resulting causal set.
        
        Implements the standard causal set sprinkling procedure:
        
        1. **Volume calculation**: V = ∫_R √|g| d^n x over region R
        2. **Poisson sampling**: N ~ Poisson(ρV) for density ρ  
        3. **Point generation**: Uniform distribution in coordinate volume
        4. **Causal matrix**: C_{ij} = 1 if x_j ≺ x_i via spacetime causality
        
        The density parameter ρ has units of inverse volume and determines
        the discretization scale. Higher densities give better continuum
        approximations but larger computational cost.

        Args:
            region (RegionDefinition): The spacetime region to sprinkle into.
            density (float): The sprinkling density ρ (points per unit volume).
            
        Returns:
            tuple: (coords, causal_matrix) where coords is an (N,d) array
                   of spacetime coordinates and causal_matrix is the (N,N)
                   causal ordering matrix.
                   
        References:
            Eq. (8) and (9): Standard Poisson sprinkling procedure with 
            N ~ Poisson(ρV) for density ρ and volume V, coordinate-independent 
            and respecting spacetime causal structure.
        """
        volume = region.get_volume(self)
        N = np.random.poisson(density * volume)

        if N == 0:
            return np.array([]).reshape(0, self.dim), np.array([]).reshape(0, 0)

        points_list = []
        while len(points_list) < N:
            candidate = region.generate_bounding_box_point(self)
            if region.contains_point(candidate, self):
                points_list.append(candidate)

        points = np.array(points_list)
        causal_matrix = np.zeros((N, N), dtype=int)
        for i in range(N):
            for j in range(N):
                if i != j and self.is_causally_related(points[j], points[i]):
                    causal_matrix[i, j] = 1

        return points, causal_matrix


def is_causally_related_minkowski(p1: np.ndarray, p2: np.ndarray) -> bool:
    """
    Checks causal precedence in flat Minkowski spacetime.
    
    Implements the Minkowski causality condition:
    p1 ≺ p2 ⟺ (Δt > 0) ∧ (Δt² - Δx² ≥ 0)
    
    where Δt = t₂ - t₁ and Δx² = |x₂ - x₁|². This corresponds to
    timelike or null separation with correct time ordering.
    
    The metric signature is (-,+,+,+) with interval:
    ds² = -dt² + dx²
    
    Args:
        p1, p2 (np.ndarray): Points in Minkowski coordinates [t, x, y, z, ...].
        
    Returns:
        bool: True if p1 causally precedes p2.
        
    References:
        Standard Minkowski causality from special relativity.
        Used as the baseline for causal set constructions in flat space.
    """
    delta = p2 - p1
    delta_t = delta[0]
    if delta_t <= 0:
        return False
    interval_sq = delta_t**2 - np.sum(delta[1:] ** 2)
    return interval_sq >= 0


class MinkowskiSpacetime(CurvedSpacetime):
    """
    Concrete implementation of 2D Minkowski spacetime for causal set theory.
    
    Implements flat spacetime with metric:
    ds² = -dt² + dx²
    
    This provides the standard testbed for causal set constructions, where
    analytic results are available for comparison. The 2D case is particularly
    tractable for Green function calculations and allows exact expressions
    for many physical quantities.
    
    **Supported Region Types:**
    - **MinkowskiDiamondRegion**: Causal diamond |t| + |x| < R
    - **MinkowskiRectangleRegion**: Rectangular region [t₁,t₂] × [x₁,x₂]  
    - **MinkowskiLightconeRegion**: Future lightcone from a spacetime point
    
    **Physical Applications:**
    - Green function calculations (exact in 2D)
    - Entanglement entropy studies
    - Spectral dimension analysis
    - Action calculations and path integrals
    
    The 2D restriction allows for exact analytical comparisons while
    retaining the essential causal structure needed for physical calculations.

    References:
        Standard 2D Minkowski spacetime. See textbooks on special relativity
        and quantum field theory in curved spacetime for the continuum theory.
    """

    def __init__(self):
        """Initialize 2D Minkowski spacetime."""
        super().__init__(dim=2)

    def is_causally_related(self, p1: np.ndarray, p2: np.ndarray) -> bool:
        """Check causality in Minkowski space."""
        return is_causally_related_minkowski(p1, p2)


class DeSitterSpacetime(CurvedSpacetime):
    """
    Concrete implementation of 2D de Sitter spacetime static patch.
    
    Implements the static patch of de Sitter space with metric:
    ds² = -(1 - r²/α²) dt² + dr²/(1 - r²/α²) + r² dΩ²
    
    In 2D coordinates (τ,χ), this becomes:
    ds² = α² dτ² / cosh²(χ) - α² dχ²
    
    where α is the de Sitter radius related to the cosmological constant
    by Λ = (d-1)(d-2)/(2α²) in d dimensions.
    
    **Coordinate Systems:**
    - **Native coordinates**: (τ, χ) covering the static patch
    - **Embedded Minkowski**: (t, x) via conformal mapping
    
    **Conformal Mapping to Minkowski:**
    The mapping to embedding Minkowski coordinates is:
    - t = α sin(τ/α) / [cosh(χ) + cos(τ/α)]
    - x = α sinh(χ) / [cosh(χ) + cos(τ/α)]
    
    **Supported Region Types:**
    - **DeSitterStaticPatchRegion**: Full static patch within horizons
    - **DeSitterDiamondRegion**: Causal diamond |τ/α ± χ| < R
    
    **Physical Applications:**
    - Cosmological causal set models
    - Horizon thermodynamics
    - de Sitter entropy calculations
    - Quantum field theory in expanding spacetime

    Attributes:
        alpha (float): The de Sitter radius parameter.
        
    References:
        Standard de Sitter geometry. See Gibbons & Hawking (1977) for
        the static patch description and thermal properties.
    """

    def __init__(self, alpha: float):
        """
        Initializes the spacetime with a de Sitter radius alpha.
        """
        super().__init__(dim=2)
        self.alpha = alpha

    def ds_to_minkowski(self, ds_coords: np.ndarray) -> np.ndarray:
        """
        Converts de Sitter (τ,χ) coordinates to embedding Minkowski (t,x).
        
        Implements the conformal mapping:
        t = α sin(τ/α) / [cosh(χ) + cos(τ/α)]  
        x = α sinh(χ) / [cosh(χ) + cos(τ/α)]
        
        This maps the static patch of de Sitter space into a region of
        2D Minkowski space, preserving the causal structure through the
        conformal factor. The mapping is bijective from the static patch
        to the interior of the diamond |t| + |x| < α.
        
        Args:
            ds_coords (np.ndarray): de Sitter coordinates [τ, χ].
            
        Returns:
            np.ndarray: Embedding Minkowski coordinates [t, x].
            
        References:
            Concept in Sec. 2 (p. 7-8): Conformal mapping from de Sitter static 
            patch to Minkowski space related to the Hawking-King-McCarthy-Malament 
            (HKMM) theorem.
        """
        tau, chi = ds_coords
        denominator = np.cosh(chi) + np.cos(tau / self.alpha)
        t = self.alpha * np.sin(tau / self.alpha) / denominator
        x = self.alpha * np.sinh(chi) / denominator
        return np.array([t, x])

    def is_causally_related(self, p1_ds: np.ndarray, p2_ds: np.ndarray) -> bool:
        """
        Checks causal precedence using conformal mapping to Minkowski space.
        
        The causal structure of de Sitter space is preserved under conformal
        transformations. This method:
        
        1. Maps both de Sitter points (τ,χ) → (t,x) via conformal embedding
        2. Applies Minkowski causality test in embedding space
        3. Returns the causal precedence relation
        
        The conformal factor cancels in the causality test since we only
        need the null cone structure, not the proper time intervals.
        
        Args:
            p1_ds, p2_ds (np.ndarray): Points in de Sitter coordinates [τ, χ].
            
        Returns:
            bool: True if p1_ds ≺ p2_ds in the de Sitter causal structure.
            
        References:
            Theorem 1 (HKMM): The causal structure determines the conformal 
            geometry. Uses conformal mapping to check causality in embedding 
            Minkowski space.
        """
        # Convert both de Sitter points to their Minkowski equivalents
        p1_mink = self.ds_to_minkowski(p1_ds)
        p2_mink = self.ds_to_minkowski(p2_ds)

        # Use the simple, reliable Minkowski checker
        return is_causally_related_minkowski(p1_mink, p2_mink)


# ============================================================================
# CAUSAL SET CLASS
# ============================================================================


class CausalSet:
    """
    Complete causal set with derived structures and physical calculations.

    This class represents a discrete spacetime as a finite partially ordered
    set, implementing the causal set approach to quantum gravity. From a set
    of spacetime coordinates, it constructs all derived causal structures
    and provides access to physical calculations.

    **Core Structures:**
    - **Causal Matrix C**: C[i,j] = 1 if x_j ≺ x_i (j in past of i)
    - **Link Matrix L**: L[i,j] = 1 if x_j → x_i is a causal link  
    - **Interval Matrix**: Number of elements in each causal interval

    **Physical Calculations Available:**
    - Green functions (massless and massive, 2D and 4D)
    - Benincasa-Dowker discrete action
    - Sorkin-Johnston vacuum and Wightman functions
    - Spacetime entanglement entropy (SSEE)
    - Spectral dimension via random walk analysis

    **Theoretical Framework:**
    The causal set represents discrete spacetime where:
    1. **Order**: Causal precedence x ≺ y encodes the light cone structure
    2. **Discreteness**: Finite cardinality provides UV regulation
    3. **Lorentz invariance**: Achieved statistically in the continuum limit

    Attributes:
        coords (np.ndarray): Spacetime coordinates (N × d array)
        N (int): Number of points in the causal set  
        causal_matrix (np.ndarray): Causal ordering matrix C (N × N)
        link_matrix (np.ndarray): Link matrix L (N × N)
        interval_sizes (np.ndarray): Interval cardinality matrix (N × N)

    References:
        Fundamental reference: Bombelli, Lee, Meyer, Sorkin (1987).
        For physical calculations: Sorkin (2007), Benincasa & Dowker (2010).
    """

    def __init__(self, coords: np.ndarray):
        """
        Initialize a causal set from coordinates.

        Args:
            coords (np.ndarray): An (N, 2) array of (t, x) coordinates.

        Raises:
            ValueError: If coords has wrong shape or contains invalid values.
        """
        if coords.ndim != 2 or (coords.size > 0 and coords.shape[1] != 2):
            msg = "Coordinates must be an (N, 2) array"
            raise ValueError(msg)

        if not np.isfinite(coords).all():
            msg = "Coordinates must contain only finite values"
            raise ValueError(msg)

        self.coords = coords.copy()  # Make a copy to avoid external modifications
        self.N = coords.shape[0]

        if self.N == 0:
            self.causal_matrix = np.array([]).reshape(0, 0)
            self.link_matrix = np.array([]).reshape(0, 0)
            self.interval_sizes = np.array([]).reshape(0, 0)
        else:
            self.causal_matrix = get_causal_matrix(coords)
            self.link_matrix = get_link_matrix(self.causal_matrix)
            self.interval_sizes = get_interval_sizes(self.causal_matrix)

    def get_Green_function_2D(self) -> np.ndarray:
        """
        Get the 2D massless scalar retarded Green function.
        
        Returns the discrete Green function:
        K₀⁽²⁾(x,x') = (1/2) C(x,x')
        
        This is exact for causal sets in 2D Minkowski spacetime and
        provides the retarded propagator for massless scalar fields.

        Returns:
            np.ndarray: The 2D Green function matrix K₀⁽²⁾.
            
        References:
            Eq. (massless2d): Exact 2D massless Green function on causal sets.
        """
        return K0_2D(self.causal_matrix)

    def get_Green_function_4D(self) -> np.ndarray:
        """
        Get the 4D massless scalar retarded Green function.
        
        Returns the discrete Green function:
        K₀⁽⁴⁾(x,x') = (1/(2π√6)) L(x,x')
        
        This approximates the 4D massless retarded Green function in the
        continuum limit ρc → ∞, where ρc is the sprinkling density.

        Returns:
            np.ndarray: The 4D Green function matrix K₀⁽⁴⁾.
            
        References:
            Eq. (eq:massless4d): 4D massless Green function approximation.
        """
        return K0_4D(self.link_matrix)

    def get_massive_Green_function(
        self, mass: float, dimension: int = 2, order: int = 2
    ) -> np.ndarray:
        """
        Get an approximation to the massive Green function.

        Args:
            mass (float): The mass of the particle.
            dimension (int): The spacetime dimension (2 or 4).
            order (int): The number of terms in the perturbative expansion.

        Returns:
            np.ndarray: The approximate massive Green function matrix.

        Raises:
            ValueError: If dimension is not 2 or 4, or mass is negative.
        """
        if mass < 0:
            msg = "Mass must be non-negative"
            raise ValueError(msg)
        if dimension not in {2, 4}:
            msg = "Dimension must be 2 or 4"
            raise ValueError(msg)

        if dimension == 2:
            G0 = self.get_Green_function_2D()
        else:
            G0 = self.get_Green_function_4D()

        return G_massive_approx(G0, mass, order)

    def calculate_BD_action(self) -> float:
        """
        Calculate the Benincasa-Dowker discrete Einstein-Hilbert action.
        
        Computes the discrete action:
        S(C) = Σₑ R(e)
        
        where R(e) is the discrete Ricci scalar at element e. This provides
        a discrete analogue of the Einstein-Hilbert action for use in
        causal set quantum gravity path integrals.

        Returns:
            float: The total Benincasa-Dowker action S(C).
            
        References:
            Benincasa & Dowker (2010): Discrete Einstein-Hilbert action.
            The action reduces to the continuum result in appropriate limits.
        """
        return calculate_BD_action(self.causal_matrix)

    def get_iDelta(self) -> np.ndarray:
        """
        Get the iΔ (Pauli-Jordan) operator for quantum field theory.
        
        Returns the discrete Pauli-Jordan function:
        iΔ = i(G_R - G_A)
        
        This operator appears in the canonical commutation relations
        [φ(x), φ(y)] = iΔ(x,y) and is fundamental for constructing
        quantum field theory on causal sets.

        Returns:
            np.ndarray: The iΔ operator matrix.
            
        References:
            Sorkin-Johnston vacuum construction using the discrete iΔ operator.
        """
        return get_iDelta_operator(self.causal_matrix)

    def get_SJ_vacuum(self) -> np.ndarray:
        """
        Get the Sorkin-Johnston vacuum Wightman function.
        
        Constructs the SJ vacuum state by minimizing the renormalized
        expectation value ⟨T(ψ²)⟩ subject to Hadamard constraints.
        
        Returns the two-point function:
        W_SJ(x,x') = ⟨0|φ(x)φ(x')|0⟩_SJ

        Returns:
            np.ndarray: The SJ vacuum Wightman function matrix.
            
        References:
            Sorkin-Johnston vacuum: Preferred vacuum state for causal sets
            that reduces to standard vacua in appropriate continuum limits.
        """
        iDelta = self.get_iDelta()
        return get_SJ_Wightman_function(iDelta)

    def calculate_entanglement_entropy(self, region_indices: list[int] | np.ndarray) -> float:
        """
        Calculate spacetime entanglement entropy (SSEE) for a subregion.
        
        Computes the entanglement entropy:
        S = -Σₖ [λₖ ln(λₖ) + (1-λₖ) ln(1-λₖ)]
        
        where λₖ are eigenvalues of the generalized eigenvalue problem
        W_A |ψₖ⟩ = λₖ iΔ_A |ψₖ⟩ restricted to the subregion A.
        
        This provides a spacetime-based generalization of entanglement
        entropy that reduces to standard results for spatial regions.

        Args:
            region_indices: Indices of points forming the subregion A.

        Returns:
            float: The entanglement entropy S of the subregion.

        Raises:
            ValueError: If region_indices contains invalid indices.
            
        References:
            Spacetime entanglement entropy following the causal set 
            generalization of the Ryu-Takayanagi proposal.
        """
        region_indices = np.asarray(region_indices)
        if region_indices.size > 0:
            if not np.all((region_indices >= 0) & (region_indices < self.N)):
                msg = "All region indices must be valid point indices"
                raise ValueError(msg)
            if len(np.unique(region_indices)) != len(region_indices):
                msg = "Region indices must be unique"
                raise ValueError(msg)

        W_SJ = self.get_SJ_vacuum()
        iDelta = self.get_iDelta()
        return calculate_ssee(W_SJ, iDelta, region_indices.tolist())

    def calculate_spectral_dimension(self, max_steps: int | None = None) -> float:
        """
        Calculate spectral dimension via discrete random walk analysis.
        
        Computes the scaling exponent:
        d_s = -2 × d(ln⟨P_σ⟩)/d(ln σ)
        
        where ⟨P_σ⟩ is the average return probability after σ random walk steps.
        The spectral dimension characterizes the effective dimensionality
        as probed by diffusion processes on the discrete geometry.

        Args:
            max_steps: Maximum random walk steps. If None, uses N//2.

        Returns:
            float: The calculated spectral dimension d_s.

        Raises:
            ValueError: If max_steps is non-positive.
            
        References:
            Spectral dimension methodology for discrete geometries.
            For continuum manifolds, d_s equals the topological dimension.
        """
        if max_steps is None:
            max_steps = max(10, self.N // 2)
        elif max_steps <= 0:
            msg = "max_steps must be positive"
            raise ValueError(msg)

        return calculate_spectral_dimension(self.link_matrix, max_steps)

    def get_physical_summary(self) -> dict:
        """
        Get a summary of key physical properties of the causal set.

        Returns:
            dict: Summary containing number of points, volume estimate,
                  action, and spectral dimension.
        """
        if self.N == 0:
            return {
                "num_points": 0,
                "volume_estimate": 0.0,
                "BD_action": 0.0,
                "spectral_dimension": np.nan,
                "num_links": 0,
            }

        # Estimate the volume from the number of points (assuming unit density)
        volume_estimate = self.N

        # Count the number of links
        num_links = np.sum(self.link_matrix)

        return {
            "num_points": self.N,
            "volume_estimate": volume_estimate,
            "BD_action": self.calculate_BD_action(),
            "spectral_dimension": self.calculate_spectral_dimension(),
            "num_links": num_links,
        }

    def __repr__(self) -> str:
        """String representation of the causal set."""
        return f"CausalSet(N={self.N}, dim=2)"

    def __len__(self) -> int:
        """Return the number of points in the causal set."""
        return self.N


# ============================================================================
# REGION DEFINITION HELPERS
# ============================================================================


def create_minkowski_diamond_region(radius: float) -> MinkowskiDiamondRegion:
    """
    Create a region definition for a Minkowski diamond.

    Args:
        radius (float): The radius of the diamond |t| + |x| < radius

    Returns:
        MinkowskiDiamondRegion: Region definition dataclass
    """
    return MinkowskiDiamondRegion(radius=radius)


def create_minkowski_rectangle_region(
    t_min: float, t_max: float, x_min: float, x_max: float
) -> MinkowskiRectangleRegion:
    """
    Create a region definition for a rectangular Minkowski region.

    Args:
        t_min, t_max (float): Time bounds
        x_min, x_max (float): Space bounds

    Returns:
        MinkowskiRectangleRegion: Region definition dataclass
    """
    return MinkowskiRectangleRegion(t_min=t_min, t_max=t_max, x_min=x_min, x_max=x_max)


def create_minkowski_lightcone_region(
    t0: float, x0: float, t_max: float
) -> MinkowskiLightconeRegion:
    """
    Create a region definition for a future lightcone in Minkowski space.

    Args:
        t0, x0 (float): Origin point of the lightcone
        t_max (float): Maximum time extent

    Returns:
        MinkowskiLightconeRegion: Region definition dataclass
    """
    return MinkowskiLightconeRegion(t0=t0, x0=x0, t_max=t_max)


def create_desitter_static_patch_region() -> DeSitterStaticPatchRegion:
    """
    Create a region definition for the full de Sitter static patch.

    Returns:
        DeSitterStaticPatchRegion: Region definition dataclass
    """
    return DeSitterStaticPatchRegion()


def create_desitter_diamond_region(radius: float) -> DeSitterDiamondRegion:
    """
    Create a region definition for a diamond region in de Sitter space.

    Args:
        radius (float): The radius parameter for the diamond

    Returns:
        DeSitterDiamondRegion: Region definition dataclass
    """
    return DeSitterDiamondRegion(radius=radius)


# ============================================================================
# USAGE EXAMPLES
# ============================================================================


def example_minkowski_regions():
    """
    Example usage of different Minkowski regions.
    """
    print("=== Minkowski Space Region Examples ===")

    # Example 1: Diamond region
    diamond_region = create_minkowski_diamond_region(radius=2.0)
    coords1 = sprinkle_minkowski_region(diamond_region, density=3.0)
    print(f"Diamond region: Generated {len(coords1)} points")

    # Example 2: Rectangular region
    rect_region = create_minkowski_rectangle_region(t_min=0, t_max=5, x_min=-2, x_max=2)
    coords2 = sprinkle_minkowski_region(rect_region, density=2.0)
    print(f"Rectangular region: Generated {len(coords2)} points")

    # Example 3: Future lightcone
    lightcone_region = create_minkowski_lightcone_region(t0=0, x0=0, t_max=3)
    coords3 = sprinkle_minkowski_region(lightcone_region, density=4.0)
    print(f"Lightcone region: Generated {len(coords3)} points")

    return coords1, coords2, coords3


def example_desitter_regions():
    """
    Example usage of different de Sitter regions.
    """
    print("=== de Sitter Space Region Examples ===")

    # Create de Sitter spacetime
    alpha = 5.0
    ds_spacetime = DeSitterSpacetime(alpha)

    # Example 1: Full static patch
    static_patch_region = create_desitter_static_patch_region()
    coords1, _causal_matrix1 = ds_spacetime.sprinkle_and_build_causet(
        static_patch_region, density=2.0
    )
    print(f"Static patch: Generated {len(coords1)} points")

    # Example 2: Diamond region in de Sitter
    diamond_region = create_desitter_diamond_region(radius=np.pi / 2)
    coords2, _causal_matrix2 = ds_spacetime.sprinkle_and_build_causet(diamond_region, density=3.0)
    print(f"de Sitter diamond: Generated {len(coords2)} points")

    return coords1, coords2


if __name__ == "__main__":
    # Run examples
    example_minkowski_regions()
    print()
    example_desitter_regions()
