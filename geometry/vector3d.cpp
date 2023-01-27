#include "vector3d.h"

template<typename F>
Vector3d<F> operator-(Vector3d<F> const &lhs, Vector3d<F> const &rhs)
{
	return Vector3d<F>{lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z};
}

template<typename F>
bool operator==(Vector3d<F> const &lhs, Vector3d<F> const &rhs)
{
	return (fabs(lhs.x-rhs.x)<1e-5) && 
	       (fabs(lhs.y-rhs.y)<1e-5) &&
	       (fabs(lhs.z-rhs.z)<1e-5);
}

template<typename F>
F operator*(Vector3d<F> const &lhs, Vector3d<F> const &rhs)
{
	return lhs.x*rhs.x+lhs.y*rhs.y+lhs.z*rhs.z;
}

template struct Vector3d<double>;
template Vector3d<double> operator- (Vector3d<double> const &, 
                                     Vector3d<double> const &);
template bool operator== (Vector3d<double> const &,
                          Vector3d<double> const &);
template double operator* (Vector3d<double> const &,
                          Vector3d<double> const &);
