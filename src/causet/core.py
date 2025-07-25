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

    Args:
        causal_matrix (np.ndarray): The N x N causal matrix.

    Returns:
        np.ndarray: The N x N link matrix L, where L[i, j] = 1
                    if j is a link to i.

    Raises:
        ValueError: If causal_matrix has wrong shape or invalid values.
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

    Args:
        causal_matrix (np.ndarray): The N x N causal matrix.

    Returns:
        np.ndarray: The N x N Green function matrix.
    """
    return 0.5 * causal_matrix


def K0_4D(link_matrix: np.ndarray) -> np.ndarray:
    """
    Calculates the 4D massless scalar Green function for a causal set.

    Args:
        link_matrix (np.ndarray): The N x N link matrix.

    Returns:
        np.ndarray: The N x N Green function matrix.
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

    Args:
        Nk_counts (dict): A dictionary with keys 0,1,2,3 for N_k(e).

    Returns:
        float: The value of R(e).
    """
    # Magic coefficients from Eq. 33
    term = 1 - Nk_counts[0] + 9 * Nk_counts[1] - 16 * Nk_counts[2] + 8 * Nk_counts[3]

    prefactor = 4 / np.sqrt(6)

    return prefactor * term


def calculate_BD_action(causal_matrix: np.ndarray) -> float:
    """
    Calculates the total Benincasa-Dowker action for a causal set.

    Args:
        causal_matrix (np.ndarray): The causal matrix for the set.

    Returns:
        float: The total action S(C).
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
    Δ = G_R - G_A. In 2D, G_R is proportional to C and G_A is proportional to C^T.
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
    Calculates the SJ Wightman function from the iΔ operator.

    Returns:
        np.ndarray: The N x N Wightman function matrix W_SJ.
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
    Calculates the Spacetime Entanglement Entropy for a subregion.

    Args:
        W_SJ (np.ndarray): The full Wightman function for the causal set.
        iDelta (np.ndarray): The full iΔ operator for the causal set.
        region_indices (list or np.ndarray): A list of indices of the points
                                             that form the subregion 'A'.

    Returns:
        float: The entanglement entropy S.
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
    Creates a transition matrix for a random walk on the causal set links.
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
    Calculates the spectral dimension by simulating a random walk.

    Args:
        link_matrix (np.ndarray): The link matrix of the causal set.
        max_steps (int): Maximum number of random walk steps to simulate.

    Returns:
        float: The calculated spectral dimension.
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
    """An abstract base class representing a general curved spacetime."""

    def __init__(self, dim: int):
        self.dim = dim

    @abstractmethod
    def is_causally_related(self, p1: np.ndarray, p2: np.ndarray) -> bool:
        """Checks for a causal relationship between two points."""

    def sprinkle_and_build_causet(self, region: RegionDefinition, density: float = 1.0):
        """Performs a Poisson sprinkling in this curved spacetime."""
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
    """Checks if point p1 causally precedes point p2 in flat Minkowski space."""
    delta = p2 - p1
    delta_t = delta[0]
    if delta_t <= 0:
        return False
    interval_sq = delta_t**2 - np.sum(delta[1:] ** 2)
    return interval_sq >= 0


class MinkowskiSpacetime(CurvedSpacetime):
    """
    A concrete implementation for 2D Minkowski spacetime.

    Supports different region types through RegionDefinition dataclasses:
    - MinkowskiDiamondRegion: Causal diamond |t| + |x| < radius
    - MinkowskiRectangleRegion: Rectangular region t_min < t < t_max, x_min < x < x_max
    - MinkowskiLightconeRegion: Future lightcone from a point
    """

    def __init__(self):
        """Initialize 2D Minkowski spacetime."""
        super().__init__(dim=2)

    def is_causally_related(self, p1: np.ndarray, p2: np.ndarray) -> bool:
        """Check causality in Minkowski space."""
        return is_causally_related_minkowski(p1, p2)


class DeSitterSpacetime(CurvedSpacetime):
    """
    A concrete implementation for a 2D de Sitter spacetime static patch.

    Coordinates:
      - Native de Sitter coords: (tau, chi)
      - Embedded Minkowski coords: (t, x)

    Supports different region types through RegionDefinition dataclasses:
    - DeSitterStaticPatchRegion: Static patch with radius bounds
    - DeSitterDiamondRegion: Causal diamond in global coordinates
    """

    def __init__(self, alpha: float):
        """
        Initializes the spacetime with a de Sitter radius alpha.
        """
        super().__init__(dim=2)
        self.alpha = alpha

    def ds_to_minkowski(self, ds_coords: np.ndarray) -> np.ndarray:
        """Converts de Sitter (tau, chi) coords to Minkowski (t, x) coords."""
        tau, chi = ds_coords
        denominator = np.cosh(chi) + np.cos(tau / self.alpha)
        t = self.alpha * np.sin(tau / self.alpha) / denominator
        x = self.alpha * np.sinh(chi) / denominator
        return np.array([t, x])

    def is_causally_related(self, p1_ds: np.ndarray, p2_ds: np.ndarray) -> bool:
        """
        Checks for causality using the conformal mapping to Minkowski space.
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
    A complete causal set with all derived properties.

    This class represents a causal set generated from a set of spacetime coordinates,
    providing access to all the derived structures and physical calculations.

    Attributes:
        coords (np.ndarray): The spacetime coordinates of the causal set points
        N (int): The number of points in the causal set
        causal_matrix (np.ndarray): The causal matrix C[i,j] = 1 if j ≺ i
        link_matrix (np.ndarray): The link matrix L[i,j] = 1 if j is a link to i
        interval_sizes (np.ndarray): Matrix of interval sizes between points
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
        Get the 2D massless scalar Green function.

        Returns:
            np.ndarray: The 2D Green function matrix.
        """
        return K0_2D(self.causal_matrix)

    def get_Green_function_4D(self) -> np.ndarray:
        """
        Get the 4D massless scalar Green function.

        Returns:
            np.ndarray: The 4D Green function matrix.
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
        Calculate the Benincasa-Dowker action.

        Returns:
            float: The total action for the causal set.
        """
        return calculate_BD_action(self.causal_matrix)

    def get_iDelta(self) -> np.ndarray:
        """
        Get the iΔ (Pauli-Jordan) operator.

        Returns:
            np.ndarray: The iΔ operator matrix.
        """
        return get_iDelta_operator(self.causal_matrix)

    def get_SJ_vacuum(self) -> np.ndarray:
        """
        Get the Sorkin-Johnston vacuum Wightman function.

        Returns:
            np.ndarray: The SJ vacuum Wightman function matrix.
        """
        iDelta = self.get_iDelta()
        return get_SJ_Wightman_function(iDelta)

    def calculate_entanglement_entropy(self, region_indices: list[int] | np.ndarray) -> float:
        """
        Calculate spacetime entanglement entropy for a subregion.

        Args:
            region_indices: Indices of points forming the subregion.

        Returns:
            float: The entanglement entropy of the subregion.

        Raises:
            ValueError: If region_indices contains invalid indices.
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
        Calculate the spectral dimension via random walk analysis.

        Args:
            max_steps: Maximum number of random walk steps. If None, uses N//2.

        Returns:
            float: The calculated spectral dimension.

        Raises:
            ValueError: If max_steps is non-positive.
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
