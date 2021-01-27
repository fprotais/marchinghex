#pragma once
#include <vector>
#include <ultimaille/geometry.h>
#include <ultimaille/volume.h>

void elliptic_smoother(UM::Hexahedra& m);
void elliptic_smoother(UM::Tetrahedra& m, std::vector<bool>& locks);






