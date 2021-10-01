#include <algorithm>
#include <cmath>
#include "Proto_BoxData.H"
#include "interpolate.h"
#include "w.h"

std::map<double (*)(double), std::tuple<int, int>> r_map = {
  {w_2, std::make_tuple(0, 1)},
  {w_3, std::make_tuple(-1, 2)},
  {w_4, std::make_tuple(-1, 2)},
  {w_6, std::make_tuple(-2, 3)},
  {L2_1, std::make_tuple(-1, 2)},
  {L2_2, std::make_tuple(-1, 2)},
  {L4_2, std::make_tuple(-2, 3)},
  {L4_4, std::make_tuple(-2, 3)}
};

vector<array<double, DIM>> compute_points_to_consider(const array<double, DIM> &original_point, const int& choice) {
  auto r_tuple = r_map[get_w(choice)];
  auto r_low = get<0>(r_tuple);
  auto r_high = get<1>(r_tuple);

  vector<array<double, DIM>> points_to_consider;

  for (int i = r_low; i <= r_high; i++) {
      for (int j = r_low; j <= r_high; j++) {
          array<double, DIM> new_1 = {original_point[0] + i, original_point[1]+j};
          points_to_consider.push_back(new_1);
      }
  }
  return points_to_consider;
}

void interpolate(BoxData<double>& grid, const double& q_k, const array<double, DIM> x_k, const double& hg, const double& hp, const int& choice)
{
  array<double, DIM> left_hand_corner = { floor(x_k[0] / hg), floor(x_k[1] / hg) };
  auto points_to_consider = compute_points_to_consider(left_hand_corner, choice);
  for (auto current_point : points_to_consider)
  {
    array<double, DIM> x_bar = { current_point[0] * hg, current_point[1] * hg };

    array<double, DIM> z = { x_bar[0] - x_k[0], x_bar[1] - x_k[1] };
    // array<double, DIM> z = { x_k[0] - x_bar[0], x_k[1] - x_bar[1] };

    double w_result = w(z, hg, choice);
    Point to_store { static_cast<int>(current_point[0]), static_cast<int>(current_point[1]) };

    // array<double, DIM> scaled_grid_point { current_point[0] * hg, current_point[1] * hg };
    if (grid.box().contains(to_store))
    {
        grid(to_store) += q_k * w_result;
    }
  }
}

// function to test whether or not the contribution of each grid point adds up to the total interpolated value
void interpolate_test(const double& q_k, const array<double, DIM> x_k, const double& hg, const double& hp, const int& choice)
{
  array<double, DIM> left_hand_corner = { floor(x_k[0] / hg), floor(x_k[1] / hg) };
  auto points_to_consider = compute_points_to_consider(left_hand_corner, choice);
  double interpolated_value = 0.;
  for (auto current_point : points_to_consider)
  {
    double contribution = 0.;
    array<double, DIM> x_bar = { current_point[0] * hg, current_point[1] * hg };

    array<double, DIM> z = { x_bar[0] - x_k[0], x_bar[1] - x_k[1] };
    // array<double, DIM> z = { x_k[0] - x_bar[0], x_k[1] - x_bar[1] };

    double w_result = w(z, hg, choice);
    Point to_store { static_cast<int>(current_point[0]), static_cast<int>(current_point[1]) };
    contribution = q_k * w_result;
    cout << "contribution: " << contribution << endl;
    interpolated_value += contribution;
  }
  cout << "full: " << interpolated_value << endl;

}
