#include "projection.h"

using namespace UM;

// Stolen from geogram
inline double distance2(const vec3& A, const vec3& B) {
    return (A - B).norm2();
}

double point_segment_squared_distance(const UM::vec3& P, const std::array<UM::vec3, 2>& AB, UM::vec3& closest_point, std::array<double, 2>& l) {
    double l2 = distance2(AB[0],AB[1]);
    double t = (P - AB[0]) * (AB[1] - AB[0]);
    if (t <= 0.0 || l2 == 0.0) {
        closest_point = AB[0];
        l[0] = 1.0;
        l[1] = 0.0;
        return distance2(P, AB[0]);
    }
    else if (t > l2) {
        closest_point = AB[1];
        l[0] = 0.0;
        l[1] = 1.0;
        return distance2(P, AB[1]);
    }
    l[1] = t / l2;
    l[0] = 1.0 - l[1];
    closest_point = l[0] * AB[0] + l[1] * AB[1];
    return distance2(P, closest_point);
}

double point_triangle_squared_distance(const UM::vec3& P, const std::array<UM::vec3, 3>& ABC, UM::vec3& closest_point, std::array<double, 3>& l) {
    vec3 diff = ABC[0] - P;
    vec3 edge0 = ABC[1] - ABC[0];
    vec3 edge1 = ABC[2] - ABC[0];
    double a00 = edge0.norm2();
    double a01 = edge0 * edge1;
    double a11 = edge1.norm2();
    double b0 = diff * edge0;
    double b1 = diff * edge1;
    double c = diff.norm2();
    double det = ::fabs(a00 * a11 - a01 * a01);
    double s = a01 * b1 - a11 * b0;
    double t = a01 * b0 - a00 * b1;
    double sqrDistance;

    // If the triangle is degenerate
    if (det < 1e-30) {
        std::array<double, 2> cur_l;
        vec3 cur_closest;
        double result;
        double cur_dist = point_segment_squared_distance(P, { ABC[0], ABC[1] }, cur_closest, cur_l);
        result = cur_dist;
        closest_point = cur_closest;
        l[0] = cur_l[0];
        l[1] = cur_l[1];
        l[2] = 0.0;
        cur_dist = point_segment_squared_distance(P, { ABC[0], ABC[2] }, cur_closest, cur_l);
        if (cur_dist < result) {
            result = cur_dist;
            closest_point = cur_closest;
            l[0] = cur_l[0];
            l[2] = cur_l[1];
            l[1] = 0.0;
        }
        cur_dist = point_segment_squared_distance(P, { ABC[1], ABC[2] }, cur_closest, cur_l);
        if (cur_dist < result) {
            result = cur_dist;
            closest_point = cur_closest;
            l[1] = cur_l[0];
            l[2] = cur_l[1];
            l[0] = 0.0;
        }
        return result;
    }

    if (s + t <= det) {
        if (s < 0.0) {
            if (t < 0.0) {   // region 4
                if (b0 < 0.0) {
                    t = 0.0;
                    if (-b0 >= a00) {
                        s = 1.0;
                        sqrDistance = a00 + 2.0 * b0 + c;
                    }
                    else {
                        s = -b0 / a00;
                        sqrDistance = b0 * s + c;
                    }
                }
                else {
                    s = 0.0;
                    if (b1 >= 0.0) {
                        t = 0.0;
                        sqrDistance = c;
                    }
                    else if (-b1 >= a11) {
                        t = 1.0;
                        sqrDistance = a11 + 2.0 * b1 + c;
                    }
                    else {
                        t = -b1 / a11;
                        sqrDistance = b1 * t + c;
                    }
                }
            }
            else {  // region 3
                s = 0.0;
                if (b1 >= 0.0) {
                    t = 0.0;
                    sqrDistance = c;
                }
                else if (-b1 >= a11) {
                    t = 1.0;
                    sqrDistance = a11 + 2.0 * b1 + c;
                }
                else {
                    t = -b1 / a11;
                    sqrDistance = b1 * t + c;
                }
            }
        }
        else if (t < 0.0) {  // region 5
            t = 0.0;
            if (b0 >= 0.0) {
                s = 0.0;
                sqrDistance = c;
            }
            else if (-b0 >= a00) {
                s = 1.0;
                sqrDistance = a00 + 2.0 * b0 + c;
            }
            else {
                s = -b0 / a00;
                sqrDistance = b0 * s + c;
            }
        }
        else {  // region 0
         // minimum at interior point
            double invDet = double(1.0) / det;
            s *= invDet;
            t *= invDet;
            sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
                t * (a01 * s + a11 * t + 2.0 * b1) + c;
        }
    }
    else {
        double tmp0, tmp1, numer, denom;

        if (s < 0.0) {   // region 2
            tmp0 = a01 + b0;
            tmp1 = a11 + b1;
            if (tmp1 > tmp0) {
                numer = tmp1 - tmp0;
                denom = a00 - 2.0 * a01 + a11;
                if (numer >= denom) {
                    s = 1.0;
                    t = 0.0;
                    sqrDistance = a00 + 2.0 * b0 + c;
                }
                else {
                    s = numer / denom;
                    t = 1.0 - s;
                    sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
                        t * (a01 * s + a11 * t + 2.0 * b1) + c;
                }
            }
            else {
                s = 0.0;
                if (tmp1 <= 0.0) {
                    t = 1.0;
                    sqrDistance = a11 + 2.0 * b1 + c;
                }
                else if (b1 >= 0.0) {
                    t = 0.0;
                    sqrDistance = c;
                }
                else {
                    t = -b1 / a11;
                    sqrDistance = b1 * t + c;
                }
            }
        }
        else if (t < 0.0) {  // region 6
            tmp0 = a01 + b1;
            tmp1 = a00 + b0;
            if (tmp1 > tmp0) {
                numer = tmp1 - tmp0;
                denom = a00 - 2.0 * a01 + a11;
                if (numer >= denom) {
                    t = 1.0;
                    s = 0.0;
                    sqrDistance = a11 + 2.0 * b1 + c;
                }
                else {
                    t = numer / denom;
                    s = 1.0 - t;
                    sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
                        t * (a01 * s + a11 * t + 2.0 * b1) + c;
                }
            }
            else {
                t = 0.0;
                if (tmp1 <= 0.0) {
                    s = 1.0;
                    sqrDistance = a00 + 2.0 * b0 + c;
                }
                else if (b0 >= 0.0) {
                    s = 0.0;
                    sqrDistance = c;
                }
                else {
                    s = -b0 / a00;
                    sqrDistance = b0 * s + c;
                }
            }
        }
        else { // region 1
            numer = a11 + b1 - a01 - b0;
            if (numer <= 0.0) {
                s = 0.0;
                t = 1.0;
                sqrDistance = a11 + 2.0 * b1 + c;
            }
            else {
                denom = a00 - 2.0 * a01 + a11;
                if (numer >= denom) {
                    s = 1.0;
                    t = 0.0;
                    sqrDistance = a00 + 2.0 * b0 + c;
                }
                else {
                    s = numer / denom;
                    t = 1.0 - s;
                    sqrDistance = s * (a00 * s + a01 * t + 2.0 * b0) +
                        t * (a01 * s + a11 * t + 2.0 * b1) + c;
                }
            }
        }
    }

    // Account for numerical round-off error.
    if (sqrDistance < 0.0) {
        sqrDistance = 0.0;
    }

    closest_point = ABC[0] + s * edge0 + t * edge1;
    l[0] = 1.0 - s - t;
    l[1] = s;
    l[2] = t;
    return sqrDistance;
}