import numpy as np
import matplotlib.pyplot as plt

# Parameters start
X = 150
Y = 150
EPSILON = 0.1
H = 0.2
C1_ALPHA = 70
C2_ALPHA = 2 * np.sqrt(C1_ALPHA)

N = 100  # Number of sensor nodes
M = 2  # Space dimensions
D = 15  # Desired distance among sensor node
K = 1.2  # Scaling factor
R = K*D  # Interaction range
DELTA_T = 0.009
A = 5
B = 5
C = np.abs(A-B)/np.sqrt(4*A*B)
ITERATION_VALUES = np.arange(0, 7, DELTA_T)
ITERATION = ITERATION_VALUES.shape[0]
SNAPSHOT_INTERVAL = 50
POSITION_X = np.zeros([N, ITERATION])
POSITION_Y = np.zeros([N, ITERATION])

# n_x = np.random.rand(N) * X
# n_y = np.random.rand(N) * X + 150
nodes = np.random.rand(N, M) * X
nodes_old = nodes
nodes_velocity_p = np.zeros([N, M])
# adjacency_matrix = np.zeros([N, N])
a_ij_matrix = np.zeros([N, N])
velocity_magnitudes = np.zeros([N, ITERATION])
connectivity = np.zeros([ITERATION, 1])
fig = plt.figure()
fig_counter = 0
c1_mt = 50
c2_mt = 2 * np.sqrt(c1_mt)
q_mt_x1 = 310
q_mt_y1 = 255
q_mt_x1_old = q_mt_x1
q_mt_y1_old = q_mt_y1
target_points = np.zeros([ITERATION, M])
center_of_mass = np.zeros([ITERATION, M])
# Parameters end

def create_adjacency_matrix():
    adjacency_matrix = np.zeros([N, N])
    for i in range(0, N):
        for j in range(0, N):
            if i != j:
                # val = nodes[i] - nodes[j]
                distance = np.linalg.norm(nodes[i] - nodes[j])
                if distance <= R:
                    adjacency_matrix[i, j] = 1
    return adjacency_matrix


def plot_deployment():
    plt.plot(nodes[:, 0], nodes[:, 1], 'ro')
    plt.show()
    # global fig_counter
    # fig_name = 'figure_' + str(fig_counter) + '.png'
    # fig_counter += 1
    # fig.savefig(fig_name, dpi=fig.dpi)


def plot_neighbors(t):
    plt.plot(target_points[0:t, 0], target_points[0:t, 1])
    # plt.plot(center_of_mass[0:t, 0], center_of_mass[0:t, 1], color='black')
    plt.plot(q_mt_x1_old, q_mt_y1_old, 'ro', color='green')
    plt.plot(nodes[:, 0], nodes[:, 1], 'ro')
    for i in range(0, N):
        for j in range(0, N):
            distance = np.linalg.norm(nodes[j] - nodes[i])
            if distance <= R:
                plt.plot([nodes[i, 0], nodes[j, 0]],
                         [nodes[i, 1], nodes[j, 1]],
                         'b-', lw=0.5)
    plt.show()


def sigma_norm(z):
    val = EPSILON*(z**2)
    val = np.sqrt(1 + val) - 1
    val = val/EPSILON
    return val


def bump_function(z):
    if 0 <= z < H:
        return 1
    elif H <= z <= 1:
        val = (z-H)/(1-H)
        val = np.cos(np.pi*val)
        val = (1+val)/2
        return val
    else:
        return 0


def sigma_1(z):
    val = 1 + z **2
    val = np.sqrt(val)
    val = z/val
    return val


def phi(z):
    val_1 = A + B
    val_2 = sigma_1(z + C)
    val_3 = A - B
    val = val_1 * val_2 + val_3
    val = val / 2
    return val


def phi_alpha(z):
    input_1 = z/sigma_norm(R)  # Sigma norm of R is R_alpha
    input_2 = z - sigma_norm(D)  # Sigma norm of D is D_alpha
    val_1 = bump_function(input_1)
    val_2 = phi(input_2)
    val = val_1 * val_2
    return val


def get_a_ij(i, j):
    val_1 = nodes[j] - nodes[i]
    norm = np.linalg.norm(val_1)
    val_2 = sigma_norm(norm)/sigma_norm(R)
    val = bump_function(val_2)
    return val


def get_n_ij(i, j):
    val_1 = nodes[j] - nodes[i]
    norm = np.linalg.norm(val_1)
    val_2 = 1 + EPSILON * norm**2
    val = val_1/np.sqrt(val_2)
    return val


def get_u_i(i, q_mt, p_mt):
    sum_1 = np.array([0.0, 0.0])
    sum_2 = np.array([0.0, 0.0])
    for j in range(0, N):
        distance = np.linalg.norm(nodes[j] - nodes[i])
        if distance <= R:
            val_1 = nodes[j] - nodes[i]
            norm = np.linalg.norm(val_1)
            phi_alpha_val = phi_alpha(sigma_norm(norm))
            val = phi_alpha_val * get_n_ij(i, j)
            sum_1 += val

            val_2 = nodes_velocity_p[j] - nodes_velocity_p[i]
            sum_2 += get_a_ij(i, j) * val_2
    val = C1_ALPHA * sum_1 + C2_ALPHA * sum_2 - c1_mt * (nodes[i] - q_mt) - \
          c2_mt * (nodes_velocity_p[i] - p_mt)
    return val


def get_positions():
    global q_mt_x1_old, q_mt_y1_old
    for t in range(0, ITERATION):
        # print(np.linalg.matrix_rank(adjacency_matrix))
        adjacency_matrix = create_adjacency_matrix()
        # print(np.linalg.matrix_rank(adjacency_matrix))
        connectivity[t] = (1 / N) * np.linalg.matrix_rank(adjacency_matrix)
        center_of_mass[t] = np.array([np.mean(nodes[:, 0]), np.mean(nodes[:, 1])])
        # print(t)
        if t == 0:
            q_mt_x1 = 310 - 160 * np.cos(ITERATION_VALUES[t])
            q_mt_y1 = 255 - 160 * np.sin(ITERATION_VALUES[t])
            q_mt = np.array([q_mt_x1, q_mt_y1])
            target_points[t] = q_mt
            plot_neighbors(t)
            for i in range(0, N):
                POSITION_X[i, t] = nodes[i, 0]
                POSITION_Y[i, t] = nodes[i, 1]
        else:
            q_mt_x1 = 310 - 160 * np.cos(ITERATION_VALUES[t])
            q_mt_y1 = 255 - 160 * np.sin(ITERATION_VALUES[t])
            q_mt = np.array([q_mt_x1, q_mt_y1])
            target_points[t] = q_mt
            q_mt_old = np.array([q_mt_x1_old, q_mt_y1_old])
            p_mt = (q_mt - q_mt_old) / DELTA_T
            q_mt_x1_old = q_mt_x1
            q_mt_y1_old = q_mt_y1
            for i in range(0, N):
                u_i = get_u_i(i, q_mt, p_mt)
                old_velocity = nodes_velocity_p[i, :]
                old_position = np.array([POSITION_X[i, t-1],
                                         POSITION_Y[i, t-1]])
                new_velocity = old_velocity + u_i * DELTA_T
                new_position = old_position + DELTA_T * new_velocity + (DELTA_T ** 2 / 2) * u_i
                # new_velocity = (new_position - old_position) / DELTA_T
                # new_velocity = old_velocity + u_i * DELTA_T
                [POSITION_X[i, t], POSITION_Y[i, t]] = new_position
                nodes_velocity_p[i, :] = new_velocity
                nodes[i, :] = new_position
                velocity_magnitudes[i, t] = np.linalg.norm(new_velocity)
                # velocity_magnitudes[i, :, t] = new_velocity
        if (t+1) % SNAPSHOT_INTERVAL == 0:
            plot_neighbors(t)


def plot_trajectory():
    for i in range(0, N):
        # arr = np.array([POSITION_X[i, :], POSITION_Y[i, :]])
        plt.plot(POSITION_X[i, :], POSITION_Y[i, :])
        # plt.show()
    plt.show()
    s = 1
    # global fig_counter
    # fig_name = 'figure_' + str(fig_counter) + '.png'
    # fig_counter += 1
    # fig.savefig(fig_name, dpi=fig.dpi)


def plot_velocity():
    for i in range(0, N):
        velocity_i = velocity_magnitudes[i, :]
        # velocity_i = 1/velocity_i
        plt.plot(velocity_i)
    plt.show()

def plot_connectivity():
    m = connectivity
    plt.plot(connectivity)
    plt.show()


def plot_center_of_mass():
    plt.plot(target_points[:, 0], target_points[:, 1])
    plt.plot(center_of_mass[:, 0], center_of_mass[:, 1])
    plt.show()


plot_deployment()
get_positions()
plot_trajectory()
plot_velocity()
plot_connectivity()
plot_center_of_mass()
