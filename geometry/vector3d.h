#pragma once
#include <cmath>

template<typename F>
struct Vector3d {
	F x, y, z;
	Vector3d<F> (F x, F y, F z): x(x), y(y), z(z) {}
	Vector3d<F> (): x(0.0), y(0.0), z(0.0) {}
};

template<typename F>
Vector3d<F> operator-(Vector3d<F> const &lhs, Vector3d<F> const &rhs);

template<typename F>
bool operator==(Vector3d<F> const &lhs, Vector3d<F> const &rhs);

template<typename F>
F operator*(Vector3d<F> const &lhs, Vector3d<F> const &rhs);
