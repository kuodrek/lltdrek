import numpy as np


def test_pool_mirroring(simple_pool):
    pool = simple_pool
    assert len(pool.pool) == 2
    assert pool.pool[0].surface_name == "main_wing"
    assert pool.pool[1].surface_name == "main_wing_mirrored"
    assert pool.pool[1].parent_wing == "main_wing"
    assert pool.pool[0].parent_wing is None
    assert np.all(pool.pool[1].collocation_points[:, 1] < 0)


def test_induced_velocities_shape(simple_pool, simple_flight_condition):
    pool = simple_pool
    for alpha in simple_flight_condition.angles_of_attack:
        for wing_i in pool.pool:
            for wing_j in pool.pool:
                v = pool.system_induced_velocities[alpha][wing_i.surface_name][wing_j.surface_name]
                assert v.shape == (wing_i.N_panels, wing_j.N_panels, 3)


def test_map_solution(simple_pool):
    pool = simple_pool
    N = pool.total_panels  # 8 = 4 original + 4 mirrored
    G_flat = np.arange(N, dtype=float)
    G_dict = pool.map_solution(G_flat)

    assert set(G_dict.keys()) == {w.surface_name for w in pool.pool}
    np.testing.assert_array_equal(G_dict["main_wing"], G_flat[:4])
    np.testing.assert_array_equal(G_dict["main_wing_mirrored"], G_flat[4:])


def test_calculate_total_velocity_shape(simple_pool, simple_flight_condition):
    pool = simple_pool
    G_dict = pool.map_solution(np.zeros(pool.total_panels))
    alpha = simple_flight_condition.angles_of_attack[0]
    total_vel = pool.calculate_total_velocity(alpha, G_dict)

    for w in pool.pool:
        assert total_vel[w.surface_name].shape == (w.N_panels, 3)

    # With G=0, total velocity equals freestream
    expected = pool.system_freestream_velocities[alpha]["main_wing"]
    np.testing.assert_allclose(total_vel["main_wing"], expected)
