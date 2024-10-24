import os
import matplotlib.pyplot as plt
from interpolation import lagrange_interpolation, generate_chebyshev_nodes, linspace, spline_interpolation

AVAILABLE_FILES = [
    ('chelm.txt', ' '),
    ('in_mountain.data', ','),
    ('stale.txt', ','),
    ('ulm_lugano.txt', ' '),
    ('tczew_starogard.txt', ' '),
    ('genoa_rapallo.txt', ' ')
]

DATA_DIR = 'data'

TEST_PARAMS = [(6, 512), (11, 512), (15, 512), (26, 512), (52, 512), (103, 512)]


def read_data(file_path, delimiter):
    with open(file_path, 'r') as file:
        return [list(map(float, line.strip().split(delimiter))) for line in file]


def generate_subplots(title, interp_x, lagrange_y, cheb_y, spline_y, selected_x, selected_y, x_data, y_data, test_params, file_name):
    fig, axs = plt.subplots(1, 3, figsize=(20, 5))

    axs[0].semilogy(interp_x, lagrange_y, label='Interpolation', color='orange')
    axs[0].semilogy(selected_x, selected_y, '.', label='Selected points', color='red')
    axs[0].semilogy(x_data, y_data, label='Input data', color='gray')
    axs[0].set_title('Lagrange Interpolation')
    axs[0].legend()

    axs[1].semilogy(interp_x, cheb_y, label='Interpolation', color='blue')
    axs[1].semilogy(selected_x, selected_y, '.', label='Selected points', color='red')
    axs[1].semilogy(x_data, y_data, label='Input data', color='gray')
    axs[1].set_title('Lagrange Interpolation with Chebyshev Nodes')
    axs[1].legend()

    axs[2].semilogy(interp_x, spline_y, label='Interpolation', color='green')
    axs[2].semilogy(selected_x, selected_y, '.', label='Selected points', color='red')
    axs[2].semilogy(x_data, y_data, label='Input data', color='gray')
    axs[2].set_title('Spline Interpolation')
    axs[2].legend()

    fig.suptitle(title)
    if not os.path.exists('plots'):
        os.makedirs('plots')

    file_dir = os.path.join('plots', os.path.splitext(os.path.basename(file_name))[0])
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)

    fig.savefig(os.path.join(file_dir, f'Interpolation_{test_params}.png'))

    plt.tight_layout()
    #plt.show()
    plt.close(fig)
    plt.close()

def select_equally_spaced_points(x_data, y_data, num_points):
    indices = linspace(0, len(x_data) - 1, num_points)
    selected_x = [x_data[int(i)] for i in indices]
    selected_y = [y_data[int(i)] for i in indices]
    return selected_x, selected_y


def select_closest_chebyshev_points(x_data, y_data, num_points):
    chebyshev_nodes = generate_chebyshev_nodes(x_data[0], x_data[-1], num_points)
    chosen_points = []
    available_points = list(zip(x_data, y_data))
    for node in chebyshev_nodes:
        closest_point = min(available_points, key=lambda point: abs(point[0] - node))
        chosen_points.append(closest_point)
        available_points.remove(closest_point)
    chosen_x = [point[0] for point in chosen_points]
    chosen_y = [point[1] for point in chosen_points]
    return chosen_x, chosen_y


def perform_interpolation_and_plot(x_data, y_data, selected_x, selected_y, file_name, num_nodes, test_params):
    interp_x = linspace(x_data[0], x_data[-1], num_nodes)
    lagrange_y = lagrange_interpolation(selected_x, selected_y, interp_x)

    cheb_x, cheb_y = select_closest_chebyshev_points(x_data, y_data, len(selected_x))
    cheb_y = lagrange_interpolation(cheb_x, cheb_y, interp_x)

    spline_y = spline_interpolation(selected_x, selected_y, interp_x)

    generate_subplots(
        f'Interpolation - {file_name} ({len(interp_x)} interpolation points, {len(selected_x)} input points)',
        interp_x, lagrange_y, cheb_y, spline_y, selected_x, selected_y, x_data, y_data, test_params, file_name
    )


def process_single_file(file_name, delimiter, test_params):
    file_path = os.path.join(DATA_DIR, file_name)
    data = read_data(file_path, delimiter)
    x_data = [point[0] for point in data]
    y_data = [point[1] for point in data]

    for points, nodes in test_params:
        selected_x, selected_y = select_equally_spaced_points(x_data, y_data, points)
        perform_interpolation_and_plot(x_data, y_data, selected_x, selected_y, file_name, nodes, points)


def choose_file():
    print(">> [INFO] Available files:")
    for file_name, _ in AVAILABLE_FILES:
        print(f"    >> {file_name}")
    print("    >> ALL (process all files)")
    while True:
        chosen_file = input(">> [INFO] Enter the name of the file you want to work on: ")
        if chosen_file.upper() == 'ALL':
            return 'ALL', None
        for file_name, delimiter in AVAILABLE_FILES:
            if file_name == chosen_file:
                return file_name, delimiter
        print(">> [ERROR] Invalid file name. Please try again.")

def main():
    file_name, delimiter = choose_file()
    if file_name == 'ALL':
        print(f'>> [INFO] Interpolating all files...')
        for file_name, delimiter in AVAILABLE_FILES:
            process_single_file(file_name, delimiter, TEST_PARAMS)
    else:
        print(f'>> [INFO] Interpolating {file_name}...')
        process_single_file(file_name, delimiter, TEST_PARAMS)
    print('>> [INFO] Done')

if __name__ == "__main__":
    main()