#include "interpolate.h"
#include "w.h"

std::map<double (*)(double), std::tuple<int, int>> r_map = {
    {w_2, std::make_tuple(0, 1)},
    {w_3, std::make_tuple(-1, 2)},
    {w_4, std::make_tuple(-1, 2)},
    {w_6, std::make_tuple(-2, 3)}};

bool point_is_present_in_vector(const array<double, DIM>& point_to_find,
                                const vector<array<double, DIM>>& vec) {
  if (find(vec.begin(), vec.end(), point_to_find) != vec.end()) {
    return true;
  }
  return false;
}

void complete_grid(vector<array<double, DIM>>& points_to_consider, const array<double, DIM>& to_add) {
  if (!point_is_present_in_vector(to_add, points_to_consider)) {
    // if the point is present, then its points to consider must already be within
    points_to_consider.push_back(to_add);
    return;
  }

  array<double, DIM> new_1 = {to_add[0] + 1, to_add[1]};
  if (!point_is_present_in_vector(new_1, points_to_consider)) {
    points_to_consider.push_back(new_1);
  }

  array<double, DIM> new_2 = {to_add[0], to_add[1] + 1};
  if (!point_is_present_in_vector(to_add, points_to_consider)) {
    points_to_consider.push_back(new_2);
  }

  array<double, DIM> new_3 = {to_add[0] + 1, to_add[1]};
  if (!point_is_present_in_vector(to_add, points_to_consider)) {
    points_to_consider.push_back(new_3);
  }
}

vector<array<double, DIM>> compute_points_to_consider(const array<double, DIM> &original_point, const int& choice) {
  auto r_tuple = r_map[get_w(choice)];
  auto r_low = get<0>(r_tuple);
  auto r_high = get<1>(r_tuple);

  vector<array<double, DIM>> points_to_consider;
  points_to_consider.push_back(original_point);

  for (int i = r_low; i <= r_high; ++i) {
    array<double, DIM> new_1 = {original_point[0] + i, original_point[1]};
    // points_to_consider.push_back(new_1);
    complete_grid(points_to_consider, new_1);
    array<double, DIM> new_2 = {original_point[0], original_point[1] + i};
    // points_to_consider.push_back(new_2);
    complete_grid(points_to_consider, new_2);
    array<double, DIM> new_3 = {original_point[0] + i, original_point[1] + i};
    // points_to_consider.push_back(new_);
    complete_grid(points_to_consider, new_3);
  }

  return points_to_consider;
}

double interpolate(const double& q_k, const array<double, DIM> x_k, const array<double, DIM> x_0, const double& h_g, const int& choice)
{
  // cout << "x_k " << x_k[0] << " , " << x_k[1] << endl;
  array<double, DIM> left_hand_corner = { floor(x_k[0] / h_g), floor(x_k[1] / h_g) };
  // array<double, DIM> left_hand_corner = { floor(x_0[0] / h_g), floor(x_0[1] / h_g) };
  cout << endl << "NEW INTERPOLATION " << endl;
  cout << "lhc: " << left_hand_corner[0] << " , " << left_hand_corner[1] << endl;
  auto points_to_consider = compute_points_to_consider(left_hand_corner, choice);
  double sum = 0;
  for (auto current_point : points_to_consider)
  {

    array<double, DIM> scaled_point = { current_point[0] * h_g, current_point[1] * h_g };
    // array<double, DIM> scaled_point = { point[0] , point[1] };

    array<double, DIM> z = { scaled_point[0] - x_k[0], scaled_point[1] - x_k[1] };
    // array<double, DIM> z = { x_k[0] - scaled_point[0], x_k[1] - scaled_point[1] };
    double w_result = w(z, h_g, choice);
    sum += q_k * w_result;
    cout << "point: " << current_point[0] << " , " << current_point[1];
    cout << "scaled point: " << scaled_point[0] << " , " << scaled_point[1] << endl;
    cout << "z: " << z[0] << " , " << z[1] << endl;
    cout << "w: " << w_result << endl;
    cout << "q_k = " << q_k << endl;
    cout << "w_result = " << w_result << endl;
  }
  cout << "sum = " << sum << endl;
  return sum;
}