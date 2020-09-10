import numpy as np
import matplotlib.pyplot as plt

# Parameters start
X = 50
Y = 50
EPSILON = 0.1
H = 0.2
C1_ALPHA = 30
C2_ALPHA = 2 * np.sqrt(C1_ALPHA)
N = 30  # Number of sensor nodes
M = 2  # Space dimensions
D = 15  # Desired distance among sensor node
K = 1.2  # Scaling factor
R = K*D  # Interaction range
DELTA_T = 0.009
A = 5
B = 5
C = np.abs(A-B)/np.sqrt(4*A*B)
ITERATION = 500
SNAPSHOT_INTERVAL = 50
POSITION_X = np.zeros([N, ITERATION])
POSITION_Y = np.zeros([N, ITERATION])


nodes = np.random.rand(N, M) * X
nodes_old = nodes
nodes_velocity_p = np.zeros([N, M])
# adjacency_matrix = np.zeros([N, N])
a_ij_matrix = np.zeros([N, N])
velocity_magnitudes = np.zeros([N, ITERATION])
connectivity = np.zeros([ITERATION, 1])
fig = plt.figure()
fig_counter = 0
q_mt = np.array([150, 150])
c1_mt = 1.1
c2_mt = 2 * np.sqrt(c1_mt)
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


def plot_neighbors():
    plt.plot(q_mt[0], q_mt[1], 'ro', color='green')
    plt.plot(nodes[:, 0], nodes[:, 1], 'ro')
    for i in range(0, N):
        for j in range(0, N):
            distance = np.linalg.norm(nodes[j] - nodes[i])
            if distance <= R:
                plt.plot([nodes[i, 0], nodes[j, 0]],
                         [nodes[i, 1], nodes[j, 1]],
                         'b-', lw=1)
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
    val_1 = nodes_old[j] - nodes_old[i]
    norm = np.linalg.norm(val_1)
    val_2 = sigma_norm(norm)/sigma_norm(R)
    val = bump_function(val_2)
    return val


def get_n_ij(i, j):
    val_1 = nodes_old[j] - nodes_old[i]
    norm = np.linalg.norm(val_1)
    val_2 = 1 + EPSILON * norm**2
    val = val_1/np.sqrt(val_2)
    return val


def get_u_i(i):
    sum_1 = np.array([0.0, 0.0])
    sum_2 = np.array([0.0, 0.0])
    for j in range(0, N):
        distance = np.linalg.norm(nodes_old[j] - nodes_old[i])
        if distance <= R:
            val_1 = nodes_old[j] - nodes_old[i]
            norm = np.linalg.norm(val_1)
            phi_alpha_val = phi_alpha(sigma_norm(norm))
            val = phi_alpha_val * get_n_ij(i, j)
            sum_1 += val

            val_2 = nodes_velocity_p[j] - nodes_velocity_p[i]
            sum_2 += get_a_ij(i, j) * val_2
    val = C1_ALPHA * sum_1 + C2_ALPHA * sum_2 - c1_mt * (nodes_old[i] - q_mt)
    return val


def get_positions():
    for t in range(0, ITERATION):
        # print(np.linalg.matrix_rank(adjacency_matrix))
        adjacency_matrix = create_adjacency_matrix()
        # print(np.linalg.matrix_rank(adjacency_matrix))
        connectivity[t] = (1 / N) * np.linalg.matrix_rank(adjacency_matrix)
        # print(t)
        if t == 0:
            plot_neighbors()
            for i in range(0, N):
                POSITION_X[i, t] = nodes[i, 0]
                POSITION_Y[i, t] = nodes[i, 1]
        else:
            if t == 3:
                o = 9
                pass
            for i in range(0, N):
                u_i = get_u_i(i)
                old_velocity = nodes_velocity_p[i, :]
                old_position = np.array([POSITION_X[i, t-1],
                                         POSITION_Y[i, t-1]])
                # new_velocity = old_velocity + u_i * DELTA_T
                new_position = old_position + DELTA_T * old_velocity + (DELTA_T ** 2 / 2) * u_i
                new_velocity = (new_position - old_position) / DELTA_T
                # new_velocity = old_velocity + u_i * DELTA_T
                [POSITION_X[i, t], POSITION_Y[i, t]] = new_position
                nodes_velocity_p[i, :] = new_velocity
                nodes[i, :] = new_position
                velocity_magnitudes[i, t] = np.linalg.norm(new_velocity)
                # velocity_magnitudes[i, :, t] = new_velocity
        if (t+1) % SNAPSHOT_INTERVAL == 0:
            plot_neighbors()


def plot_trajectory():
    for i in range(0, N):
        plt.plot(POSITION_X[i, :], POSITION_Y[i, :])

    plt.show()



def plot_velocity():
    for i in range(0, N):
        velocity_i = velocity_magnitudes[i, :]
        plt.plot(velocity_i)
    plt.show()



def plot_connectivity():
    plt.plot(connectivity)
    plt.show()


plot_deployment()
get_positions()
plot_trajectory()
plot_velocity()
plot_connectivity()

