"""
Test suite for the causal set theory core module.

Tests the dataclass-based region system, sprinkling algorithms,
causal matrix computation, and all related functionality.
"""

import warnings

import numpy as np
import pytest

from src.causet.core import (
    calculate_BD_action,
    # Causal set functionality
    CausalSet,
    create_desitter_diamond_region,
    create_desitter_static_patch_region,
    # Helper functions
    create_minkowski_diamond_region,
    create_minkowski_lightcone_region,
    create_minkowski_rectangle_region,
    DeSitterDiamondRegion,
    DeSitterSpacetime,
    DeSitterStaticPatchRegion,
    G_massive_approx,
    get_causal_matrix,
    get_iDelta_operator,
    get_interval_sizes,
    get_link_matrix,
    get_SJ_Wightman_function,
    is_causally_related_minkowski,
    # Green functions and actions (using correct names)
    K0_2D,
    K0_4D,
    # Region dataclasses
    MinkowskiDiamondRegion,
    MinkowskiLightconeRegion,
    MinkowskiRectangleRegion,
    # Spacetime classes
    MinkowskiSpacetime,
    sprinkle_minkowski_diamond,
    # Sprinkling functions
    sprinkle_minkowski_region,
)


class TestRegionDataclasses:
    """Test the region definition dataclasses."""

    def test_minkowski_diamond_region(self):
        """Test MinkowskiDiamondRegion functionality."""
        radius = 2.0
        diamond = MinkowskiDiamondRegion(radius=radius)
        spacetime = MinkowskiSpacetime()  # Need spacetime for methods

        # Test volume calculation
        expected_volume = 2 * radius**2
        assert abs(diamond.get_volume(spacetime) - expected_volume) < 1e-10

        # Test point generation
        point = diamond.generate_bounding_box_point(spacetime)
        assert len(point) == 2
        assert -radius <= point[0] <= radius
        assert -radius <= point[1] <= radius

        # Test point containment
        # Point inside diamond
        inside_point = np.array([0.5, 0.5])
        assert diamond.contains_point(inside_point, spacetime)

        # Point outside diamond
        outside_point = np.array([1.5, 1.5])
        assert not diamond.contains_point(outside_point, spacetime)

        # Point on boundary (should be outside due to strict inequality)
        boundary_point = np.array([1.0, 1.0])
        assert not diamond.contains_point(boundary_point, spacetime)

        # Test compatibility
        assert diamond.validate_spacetime_compatibility(spacetime)

    def test_minkowski_diamond_region_validation(self):
        """Test MinkowskiDiamondRegion input validation."""
        # Negative radius
        with pytest.raises(ValueError, match="Radius must be positive"):
            MinkowskiDiamondRegion(radius=-1.0)

        # Test wrong spacetime dimension
        diamond = MinkowskiDiamondRegion(radius=1.0)

        class Mock3DSpacetime:
            dim = 3

        with pytest.raises(ValueError, match="only supports 2D spacetime"):
            diamond.get_volume(Mock3DSpacetime())

        # Test wrong point dimension
        spacetime = MinkowskiSpacetime()
        with pytest.raises(ValueError, match="Point must be 2D"):
            diamond.contains_point(np.array([1, 2, 3]), spacetime)

    def test_minkowski_rectangle_region(self):
        """Test MinkowskiRectangleRegion functionality."""
        t_min, t_max = 0.0, 2.0
        x_min, x_max = -1.0, 1.0
        rect = MinkowskiRectangleRegion(t_min=t_min, t_max=t_max, x_min=x_min, x_max=x_max)
        spacetime = MinkowskiSpacetime()

        # Test volume calculation
        expected_volume = (t_max - t_min) * (x_max - x_min)
        assert abs(rect.get_volume(spacetime) - expected_volume) < 1e-10

        # Test point generation
        point = rect.generate_bounding_box_point(spacetime)
        assert len(point) == 2
        assert t_min <= point[0] <= t_max
        assert x_min <= point[1] <= x_max

        # Test point containment
        inside_point = np.array([1.0, 0.0])
        assert rect.contains_point(inside_point, spacetime)

        outside_point = np.array([3.0, 0.0])
        assert not rect.contains_point(outside_point, spacetime)

        # Test compatibility
        assert rect.validate_spacetime_compatibility(spacetime)

    def test_minkowski_rectangle_region_validation(self):
        """Test MinkowskiRectangleRegion input validation."""
        # Invalid time bounds
        with pytest.raises(ValueError, match="t_max must be greater than t_min"):
            MinkowskiRectangleRegion(t_min=2.0, t_max=1.0, x_min=-1.0, x_max=1.0)

        # Invalid space bounds
        with pytest.raises(ValueError, match="x_max must be greater than x_min"):
            MinkowskiRectangleRegion(t_min=0.0, t_max=1.0, x_min=1.0, x_max=-1.0)

    def test_minkowski_lightcone_region(self):
        """Test MinkowskiLightconeRegion functionality."""
        t0, x0, t_max = 0.0, 0.0, 2.0
        lightcone = MinkowskiLightconeRegion(t0=t0, x0=x0, t_max=t_max)
        spacetime = MinkowskiSpacetime()

        # Test volume calculation
        height = t_max - t0
        expected_volume = height**2
        assert abs(lightcone.get_volume(spacetime) - expected_volume) < 1e-10

        # Test point generation
        point = lightcone.generate_bounding_box_point(spacetime)
        assert len(point) == 2

        # Test point containment
        # Point inside lightcone
        inside_point = np.array([1.0, 0.5])
        assert lightcone.contains_point(inside_point, spacetime)

        # Point outside lightcone (too early)
        outside_point1 = np.array([-0.5, 0.0])
        assert not lightcone.contains_point(outside_point1, spacetime)

        # Point outside lightcone (too far spatially)
        outside_point2 = np.array([1.0, 1.5])
        assert not lightcone.contains_point(outside_point2, spacetime)

        # Test compatibility
        assert lightcone.validate_spacetime_compatibility(spacetime)

    def test_minkowski_lightcone_region_validation(self):
        """Test MinkowskiLightconeRegion input validation."""
        # Invalid time bounds
        with pytest.raises(ValueError, match="t_max must be greater than t0"):
            MinkowskiLightconeRegion(t0=2.0, x0=0.0, t_max=1.0)

    def test_desitter_static_patch_region(self):
        """Test DeSitterStaticPatchRegion functionality."""
        patch = DeSitterStaticPatchRegion()
        alpha = 5.0
        spacetime = DeSitterSpacetime(alpha)

        # Test volume calculation
        expected_volume = alpha**2 * np.pi**2
        assert abs(patch.get_volume(spacetime) - expected_volume) < 1e-10

        # Test point generation
        point = patch.generate_bounding_box_point(spacetime)
        assert len(point) == 2

        # Test that generated points are always in region (by design)
        assert patch.contains_point(point, spacetime)

        # Test compatibility
        assert patch.validate_spacetime_compatibility(spacetime)

    def test_desitter_diamond_region(self):
        """Test DeSitterDiamondRegion functionality."""
        radius = np.pi / 2
        diamond = DeSitterDiamondRegion(radius=radius)
        alpha = 5.0
        spacetime = DeSitterSpacetime(alpha)

        # Test volume calculation
        expected_volume = alpha**2 * radius**2
        assert abs(diamond.get_volume(spacetime) - expected_volume) < 1e-10

        # Test point generation
        point = diamond.generate_bounding_box_point(spacetime)
        assert len(point) == 2

        # Test compatibility
        assert diamond.validate_spacetime_compatibility(spacetime)

    def test_desitter_diamond_region_validation(self):
        """Test DeSitterDiamondRegion input validation."""
        # Negative radius
        with pytest.raises(ValueError, match="Radius must be positive"):
            DeSitterDiamondRegion(radius=-1.0)

        # Large radius warning (should not raise but gives warning)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            DeSitterDiamondRegion(radius=2 * np.pi)
            assert len(w) == 1
            assert "may exceed static patch boundaries" in str(w[0].message)


class TestHelperFunctions:
    """Test the helper functions that create region instances."""

    def test_create_minkowski_diamond_region(self):
        """Test create_minkowski_diamond_region helper."""
        radius = 3.0
        diamond = create_minkowski_diamond_region(radius)

        assert isinstance(diamond, MinkowskiDiamondRegion)
        assert diamond.radius == radius

    def test_create_minkowski_rectangle_region(self):
        """Test create_minkowski_rectangle_region helper."""
        t_min, t_max = 1.0, 3.0
        x_min, x_max = -2.0, 2.0
        rect = create_minkowski_rectangle_region(t_min, t_max, x_min, x_max)

        assert isinstance(rect, MinkowskiRectangleRegion)
        assert rect.t_min == t_min
        assert rect.t_max == t_max
        assert rect.x_min == x_min
        assert rect.x_max == x_max

    def test_create_minkowski_lightcone_region(self):
        """Test create_minkowski_lightcone_region helper."""
        t0, x0, t_max = 1.0, 0.5, 4.0
        lightcone = create_minkowski_lightcone_region(t0, x0, t_max)

        assert isinstance(lightcone, MinkowskiLightconeRegion)
        assert lightcone.t0 == t0
        assert lightcone.x0 == x0
        assert lightcone.t_max == t_max

    def test_create_desitter_static_patch_region(self):
        """Test create_desitter_static_patch_region helper."""
        patch = create_desitter_static_patch_region()

        assert isinstance(patch, DeSitterStaticPatchRegion)

    def test_create_desitter_diamond_region(self):
        """Test create_desitter_diamond_region helper."""
        radius = np.pi / 3
        diamond = create_desitter_diamond_region(radius)

        assert isinstance(diamond, DeSitterDiamondRegion)
        assert diamond.radius == radius


class TestSprinklingFunctions:
    """Test the point sprinkling algorithms."""

    def test_sprinkle_minkowski_region_diamond(self):
        """Test sprinkling in a Minkowski diamond region."""
        radius = 2.0
        density = 5.0
        diamond = create_minkowski_diamond_region(radius)
        spacetime = MinkowskiSpacetime()

        coords = sprinkle_minkowski_region(diamond, density)

        # Check that we got a reasonable number of points
        expected_count = density * diamond.get_volume(spacetime)
        assert len(coords) > 0
        assert abs(len(coords) - expected_count) < 3 * np.sqrt(expected_count)  # 3-sigma

        # Check that all points are in the region
        for point in coords:
            assert diamond.contains_point(point, spacetime)

        # Check coordinate structure
        assert coords.shape[1] == 2  # 2D coordinates

    def test_sprinkle_minkowski_region_rectangle(self):
        """Test sprinkling in a rectangular region."""
        rect = create_minkowski_rectangle_region(0, 2, -1, 1)
        density = 3.0
        spacetime = MinkowskiSpacetime()

        coords = sprinkle_minkowski_region(rect, density)

        # Check that all points are in the region
        for point in coords:
            assert rect.contains_point(point, spacetime)

    def test_sprinkle_minkowski_region_lightcone(self):
        """Test sprinkling in a lightcone region."""
        lightcone = create_minkowski_lightcone_region(0, 0, 2)
        density = 4.0
        spacetime = MinkowskiSpacetime()

        coords = sprinkle_minkowski_region(lightcone, density)

        # Check that all points are in the region
        for point in coords:
            assert lightcone.contains_point(point, spacetime)

    def test_sprinkle_minkowski_diamond_compatibility(self):
        """Test that sprinkle_minkowski_diamond still works."""
        radius = 1.5
        density = 3.0
        spacetime = MinkowskiSpacetime()

        coords = sprinkle_minkowski_diamond(radius, density)

        # Check that all points are in the diamond
        diamond = MinkowskiDiamondRegion(radius=radius)
        for point in coords:
            assert diamond.contains_point(point, spacetime)

    def test_sprinkle_minkowski_diamond_validation(self):
        """Test validation in sprinkle_minkowski_diamond."""
        # Negative radius
        with pytest.raises(ValueError, match="Radius must be positive"):
            sprinkle_minkowski_diamond(radius=-1.0)

        # Zero density
        with pytest.raises(ValueError, match="Density must be positive"):
            sprinkle_minkowski_diamond(radius=1.0, density=0.0)

    def test_sprinkle_minkowski_region_validation(self):
        """Test validation in sprinkle_minkowski_region."""
        diamond = create_minkowski_diamond_region(1.0)

        # Negative density
        with pytest.raises(ValueError, match="Density must be positive"):
            sprinkle_minkowski_region(diamond, density=-1.0)

        # Incompatible region (simulate by using a de Sitter region)
        ds_region = create_desitter_static_patch_region()
        with pytest.raises(ValueError, match="is not compatible with MinkowskiSpacetime"):
            sprinkle_minkowski_region(ds_region, density=1.0)


class TestSpacetimeClasses:
    """Test the spacetime implementation classes."""

    def test_minkowski_spacetime_creation(self):
        """Test MinkowskiSpacetime initialization."""
        spacetime = MinkowskiSpacetime()
        assert spacetime.dim == 2

    def test_minkowski_spacetime_causality(self):
        """Test causality checking in Minkowski spacetime."""
        spacetime = MinkowskiSpacetime()

        # Test timelike separation (causal)
        p1 = np.array([0.0, 0.0])
        p2 = np.array([2.0, 1.0])
        assert spacetime.is_causally_related(p1, p2)

        # Test spacelike separation (not causal)
        p3 = np.array([0.0, 0.0])
        p4 = np.array([1.0, 2.0])
        assert not spacetime.is_causally_related(p3, p4)

    def test_desitter_spacetime_creation(self):
        """Test DeSitterSpacetime initialization."""
        alpha = 3.0
        spacetime = DeSitterSpacetime(alpha)
        assert spacetime.dim == 2
        assert spacetime.alpha == alpha

    def test_desitter_coordinate_conversion(self):
        """Test de Sitter to Minkowski coordinate conversion."""
        alpha = 2.0
        spacetime = DeSitterSpacetime(alpha)

        # Test coordinate conversion
        ds_coords = np.array([0.0, 0.0])
        mink_coords = spacetime.ds_to_minkowski(ds_coords)
        assert len(mink_coords) == 2
        assert np.isfinite(mink_coords).all()


class TestCausalSetFunctionality:
    """Test causal set creation and analysis."""

    def test_causal_set_creation(self):
        """Test CausalSet creation from coordinates."""
        # Create simple test coordinates
        coords = np.array([[0.0, 0.0], [1.0, 0.5], [2.0, 0.0]])

        causet = CausalSet(coords)
        assert causet.N == 3
        assert np.array_equal(causet.coords, coords)
        assert len(causet) == 3
        assert str(causet) == "CausalSet(N=3, dim=2)"

    def test_causal_set_empty(self):
        """Test CausalSet creation with empty coordinates."""
        coords = np.array([]).reshape(0, 2)
        causet = CausalSet(coords)
        assert causet.N == 0
        assert len(causet) == 0
        assert causet.causal_matrix.shape == (0, 0)

    def test_causal_set_validation(self):
        """Test CausalSet input validation."""
        # Wrong shape
        with pytest.raises(ValueError, match="Coordinates must be an \\(N, 2\\) array"):
            CausalSet(np.array([1, 2, 3]))

        # Wrong number of dimensions
        with pytest.raises(ValueError, match="Coordinates must be an \\(N, 2\\) array"):
            CausalSet(np.array([[1], [2]]))

        # Non-finite values
        with pytest.raises(ValueError, match="Coordinates must contain only finite values"):
            CausalSet(np.array([[0, 0], [np.inf, 0]]))

    def test_causal_set_methods(self):
        """Test various CausalSet methods."""
        coords = np.array([[0.0, 0.0], [1.0, 0.5], [2.0, 0.0]])

        causet = CausalSet(coords)

        # Test Green functions
        G2D = causet.get_Green_function_2D()
        assert G2D.shape == (3, 3)

        G4D = causet.get_Green_function_4D()
        assert G4D.shape == (3, 3)

        # Test massive Green function
        Gm = causet.get_massive_Green_function(mass=0.1, dimension=2)
        assert Gm.shape == (3, 3)

        # Test action
        action = causet.calculate_BD_action()
        assert np.isfinite(action)

        # Test spectral dimension
        spec_dim = causet.calculate_spectral_dimension(max_steps=5)
        assert np.isfinite(spec_dim) or np.isnan(spec_dim)

        # Test physical summary
        summary = causet.get_physical_summary()
        assert "num_points" in summary
        assert summary["num_points"] == 3

    def test_causal_set_massive_green_validation(self):
        """Test validation in massive Green function."""
        coords = np.array([[0, 0], [1, 0]])
        causet = CausalSet(coords)

        # Negative mass
        with pytest.raises(ValueError, match="Mass must be non-negative"):
            causet.get_massive_Green_function(mass=-1.0)

        # Invalid dimension
        with pytest.raises(ValueError, match="Dimension must be 2 or 4"):
            causet.get_massive_Green_function(mass=1.0, dimension=3)

    def test_causal_set_entanglement_entropy_validation(self):
        """Test validation in entanglement entropy calculation."""
        coords = np.array([[0, 0], [1, 0], [2, 0]])
        causet = CausalSet(coords)

        # Invalid indices
        with pytest.raises(ValueError, match="All region indices must be valid"):
            causet.calculate_entanglement_entropy([0, 5])

        # Duplicate indices
        with pytest.raises(ValueError, match="Region indices must be unique"):
            causet.calculate_entanglement_entropy([0, 0, 1])

        # Valid calculation
        entropy = causet.calculate_entanglement_entropy([0, 1])
        assert np.isfinite(entropy)

    def test_causal_set_spectral_dimension_validation(self):
        """Test validation in spectral dimension calculation."""
        coords = np.array([[0, 0], [1, 0]])
        causet = CausalSet(coords)

        # Invalid max_steps
        with pytest.raises(ValueError, match="max_steps must be positive"):
            causet.calculate_spectral_dimension(max_steps=0)

    def test_is_causally_related_minkowski(self):
        """Test the basic causality function."""
        # Timelike separation
        p1 = np.array([0.0, 0.0])
        p2 = np.array([2.0, 1.0])
        assert is_causally_related_minkowski(p1, p2)

        # Spacelike separation
        p3 = np.array([0.0, 0.0])
        p4 = np.array([1.0, 2.0])
        assert not is_causally_related_minkowski(p3, p4)

        # Lightlike separation (boundary case)
        p5 = np.array([0.0, 0.0])
        p6 = np.array([1.0, 1.0])
        # Should be True (lightlike is causal)
        assert is_causally_related_minkowski(p5, p6)

        # Backwards in time
        p7 = np.array([1.0, 0.0])
        p8 = np.array([0.0, 0.0])
        assert not is_causally_related_minkowski(p7, p8)


class TestGreenFunctionsAndActions:
    """Test Green functions and action calculations."""

    def test_K0_2D_green_function(self):
        """Test 2D Green function calculation."""
        # Create a simple causal matrix
        causal_matrix = np.array([[0, 1, 1], [0, 0, 1], [0, 0, 0]])

        G = K0_2D(causal_matrix)
        assert G.shape == (3, 3)
        assert np.isfinite(G).all()
        assert np.array_equal(G, 0.5 * causal_matrix)

    def test_K0_4D_green_function(self):
        """Test 4D Green function calculation."""
        # Create a simple link matrix
        link_matrix = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]])

        G = K0_4D(link_matrix)
        assert G.shape == (3, 3)
        assert np.isfinite(G).all()

        expected_prefactor = 1 / (2 * np.pi * np.sqrt(6))
        assert np.allclose(G, expected_prefactor * link_matrix)

    def test_massive_green_function(self):
        """Test massive Green function approximation."""
        # Create a simple Green function matrix
        G0 = np.array([[0.0, 0.5, 0.5], [0.0, 0.0, 0.5], [0.0, 0.0, 0.0]])

        mass = 0.1
        order = 2
        Gm = G_massive_approx(G0, mass, order)

        assert Gm.shape == G0.shape
        assert np.isfinite(Gm).all()

        # For small mass, should be close to massless case
        assert np.allclose(Gm, G0, atol=0.1)

    def test_causal_matrix_validation(self):
        """Test causal matrix input validation."""
        # Empty coordinates
        coords = np.array([]).reshape(0, 2)
        C = get_causal_matrix(coords)
        assert C.shape == (0, 0)

        # Wrong shape
        with pytest.raises(ValueError, match="Coordinates must be an \\(N, 2\\) array"):
            get_causal_matrix(np.array([1, 2, 3]))

        # Non-finite values
        with pytest.raises(ValueError, match="Coordinates must contain only finite values"):
            get_causal_matrix(np.array([[0, 0], [np.inf, 0]]))

    def test_link_matrix_validation(self):
        """Test link matrix input validation."""
        # Non-square matrix
        with pytest.raises(ValueError, match="Causal matrix must be square"):
            get_link_matrix(np.array([[0, 1], [1, 0], [0, 0]]))

        # Invalid values
        with pytest.raises(ValueError, match="Causal matrix must contain only 0s and 1s"):
            get_link_matrix(np.array([[0, 2], [0, 0]]))

    def test_interval_sizes_validation(self):
        """Test interval sizes input validation."""
        # Non-square matrix
        with pytest.raises(ValueError, match="Causal matrix must be square"):
            get_interval_sizes(np.array([[0, 1], [1, 0], [0, 0]]))

        # Invalid values
        with pytest.raises(ValueError, match="Causal matrix must contain only 0s and 1s"):
            get_interval_sizes(np.array([[0, 2], [0, 0]]))

    def test_sj_wightman_function_validation(self):
        """Test SJ Wightman function input validation."""
        # Non-square matrix
        with pytest.raises(ValueError, match="iDelta operator must be square"):
            get_SJ_Wightman_function(np.array([[0, 1], [1, 0], [0, 0]]))

    def test_calculate_BD_action(self):
        """Test Benincasa-Dowker action calculation."""
        # Create a simple causal matrix
        causal_matrix = np.array([[0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1], [0, 0, 0, 0]])

        # Calculate action
        action = calculate_BD_action(causal_matrix)
        assert np.isfinite(action)
        assert isinstance(action, (float, np.floating))

    def test_sj_wightman_function(self):
        """Test SJ Wightman function calculation."""
        # Create test causal matrix
        causal_matrix = np.array([[0, 1, 1], [0, 0, 1], [0, 0, 0]])

        # Get iDelta operator first
        iDelta = get_iDelta_operator(causal_matrix)

        # Calculate Wightman function
        W_SJ = get_SJ_Wightman_function(iDelta)
        assert W_SJ.shape == causal_matrix.shape
        assert np.isfinite(W_SJ).all()


class TestIntegrationExamples:
    """Integration tests using realistic examples."""

    def test_full_workflow_minkowski(self):
        """Test complete workflow in Minkowski spacetime."""
        # Create spacetime and region
        MinkowskiSpacetime()
        diamond = create_minkowski_diamond_region(radius=2.0)

        # Sprinkle points
        coords = sprinkle_minkowski_region(diamond, density=5.0)
        assert len(coords) > 0

        # Create causal set
        causet = CausalSet(coords)
        assert causet.N == len(coords)

        # Test some physics calculations
        if causet.N >= 2:
            # Calculate causal matrix
            causal_matrix = get_causal_matrix(coords)
            assert causal_matrix.shape == (causet.N, causet.N)

            # Calculate Green function matrix
            G_matrix = K0_2D(causal_matrix)
            assert G_matrix.shape == (causet.N, causet.N)
            assert np.isfinite(G_matrix).all()

            # Calculate action
            action = calculate_BD_action(causal_matrix)
            assert np.isfinite(action)

    def test_full_workflow_desitter(self):
        """Test complete workflow in de Sitter spacetime."""
        # Create spacetime and region
        alpha = 3.0
        spacetime = DeSitterSpacetime(alpha)
        patch = create_desitter_static_patch_region()

        # Use the abstract interface
        coords, causal_matrix = spacetime.sprinkle_and_build_causet(patch, density=2.0)
        assert len(coords) > 0
        assert causal_matrix.shape == (len(coords), len(coords))

        # Create causal set
        causet = CausalSet(coords)
        assert causet.N == len(coords)

    def test_multiple_region_types(self):
        """Test that all region types work correctly."""
        spacetime = MinkowskiSpacetime()
        regions = [
            create_minkowski_diamond_region(radius=1.5),
            create_minkowski_rectangle_region(0, 2, -1, 1),
            create_minkowski_lightcone_region(0, 0, 1.5),
        ]

        for region in regions:
            coords = sprinkle_minkowski_region(region, density=3.0)
            assert len(coords) > 0

            # Verify all points are in region
            for point in coords:
                assert region.contains_point(point, spacetime)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
