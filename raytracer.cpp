#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
#include <initializer_list>
#include <string>
#include <algorithm>
#include <tuple>
#include <list>
#include <random>
#include <chrono>

static std::uniform_real_distribution<double> uniform(0, 1);
static std::default_random_engine engine(10);

class Vector
{
public:
	explicit Vector(double x = 0, double y = 0, double z = 0)
	{
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	double norm2() const
	{
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() const
	{
		return sqrt(norm2());
	}
	void normalize()
	{
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
	double operator[](int i) const { return data[i]; };
	double &operator[](int i) { return data[i]; };
	double data[3];
};

Vector operator+(const Vector &a, const Vector &b)
{
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b)
{
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator-(const Vector &a)
{
	return Vector(-a[0], -a[1], -a[2]);
}
Vector operator*(const double a, const Vector &b)
{
	return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &a, const double b)
{
	return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator/(const Vector &a, const double b)
{
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}

double dot(const Vector &a, const Vector &b)
{

	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector &a, const Vector &b)
{
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
Vector operator*(const Vector a, const Vector &b)
{
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

class Ray
{
public:
	explicit Ray(Vector v, Vector u)
	{
		origin = v;
		direction = u;
	}

	Vector origin;
	Vector direction;
};

// done in class
Vector random_cos(Vector &N)
{
	double r1 = uniform(engine);
	double r2 = uniform(engine);
	double s = sqrt(1 - r2);
	double x = cos(2 * M_PI * r1) * s;
	double y = sin(2 * M_PI * r1) * s;
	double z = sqrt(r2);

	Vector T1;
	if (abs(N[1]) < abs(N[0]) && abs(N[1]) < abs(N[2]))
	{
		T1 = Vector(N[2], 0, -N[0]);
	}
	else if (abs(N[0]) < abs(N[1]) && abs(N[0]) < abs(N[2]))
	{
		T1 = Vector(0, N[2], -N[1]);
	}
	else
	{
		T1 = Vector(N[1], -N[0], 0);
	}
	T1.normalize();
	Vector T2 = cross(N, T1);

	Vector V = x * T1 + y * T2 + z * N;
	return V;
}

class Object
{
public:
	Object(const Vector &albedo, bool mirror = false, bool transparent = false)
	{
		this->albedo = albedo;
		this->mirror = mirror;
		this->transparent = transparent;
	}
	virtual bool intersection(const Ray &r, Vector &P, Vector &N, double &t) = 0;

	Vector albedo;
	bool mirror, transparent;
};

class BoundingBox
{
public:
	BoundingBox(){};
	bool rayBoxIntersection(const Ray &ray)
	{
		Vector t0(1E9, 1E9, 1E9);
		Vector t1(-1E9, -1E9, -1E9);
		double pt1, pt2;
		for (int i = 0; i < 3; ++i)
		{
			Vector N((i == 0), (i == 1), (i == 2));
			double uDotN = dot(ray.direction, N);
			pt1 = dot(m - ray.origin, N) / uDotN;
			pt2 = dot(M - ray.origin, N) / uDotN;
			t0[i] = std::min(pt1, pt2);
			t1[i] = std::max(pt1, pt2);
		}
		double t1_min = std::min(std::min(t1[0], t1[1]), t1[2]);
		double t0_min = std::max(std::max(t0[0], t0[1]), t0[2]);
		return (t1_min > t0_min > 0);
	}
	Vector m, M;
};

class Sphere : public Object
{
public:
	explicit Sphere(Vector center, double radius, const Vector a, bool m = false, bool t = false) : center(center), radius(radius), Object(a, m, t) {}

	bool intersection(const Ray &r, Vector &P, Vector &N, double &t) override
	{
		Vector ominusc = r.origin - (this->center);
		double delta = dot(r.direction, ominusc) * dot(r.direction, ominusc) - ominusc.norm2() + this->radius * this->radius;
		if (delta >= 0)
		{
			double t1 = -dot(r.direction, ominusc) - sqrt(delta);
			double t2 = -dot(r.direction, ominusc) + sqrt(delta);
			if (t2 < 0)
				return false;
			if (t1 > 0)
			{
				t = t1;
			}
			else
			{
				t = t2;
			}
			P = r.origin + t * r.direction;
			N = P - center;
			N.normalize();
			return true;
		}

		return false;
	}

	Vector center;
	double radius;
};

double sqr(double a)
{
	return a * a;
}

class Scene
{
public:
	Scene(const double I, Vector scene_light, std::vector<Object *> objects = {}) : I(I), scene_light(scene_light), objects(objects){};

	void addToScene(Object *sphere)
	{
		objects.push_back(sphere);
	}

	bool sceneIntersection(const Ray &ray, double &t, Vector &P, Vector &N, int &id)
	{
		double t_object;
		Vector P_object, N_object;
		bool inter = false;
		t = 1E9;
		for (int i = 0; i < objects.size(); i++)
		{
			bool object_inter = objects[i]->intersection(ray, P_object, N_object, t_object);
			inter = inter || object_inter;
			if (object_inter && (t_object < t))
			{
				t = t_object;
				P = P_object;
				N = N_object;
				id = i;
			}
		}
		return inter;
	}

	Vector getColor(const Ray &ray, int bounce)
	{
		Vector color_pixel(0, 0, 0);
		if (bounce <= 0)
			return color_pixel;
		Vector P, N;
		int indexElement;
		double t;

		if (sceneIntersection(ray, t, P, N, indexElement))
		{

			if (objects[indexElement]->mirror)
			{
				Vector R = ray.direction - 2 * dot(ray.direction, N) * N;
				Ray reflectionRay(P + 0.001 * N, R);
				return getColor(reflectionRay, bounce - 1);
			}

			if (objects[indexElement]->transparent)
			{
				Vector R = ray.direction - 2 * dot(ray.direction, N) * N;
				Ray reflectionRay(P + 0.001 * N, R);
				double n1 = 1, n2 = 1.4;
				Vector NforTransp = N;
				if (dot(ray.direction, N) > 0)
				{
					std::swap(n1, n2);
					NforTransp = -NforTransp;
				}
				Vector Tt = n1 / n2 * (ray.direction - dot(ray.direction, N) * N);
				double radic = 1 - ((n1 / n2) * (n1 / n2)) * (1 - (dot(ray.direction, N) * dot(ray.direction, N)));
				if (radic < 0)
					return getColor(reflectionRay, bounce - 1);
				Vector Tn = -sqrt(radic) * NforTransp;
				Vector T = Tt + Tn;
				Ray refractionRay(P + 0.001 * T, T);
				return getColor(refractionRay, bounce - 1);
			}

			Vector light = scene_light - P;
			double distlight = light.norm2();
			light.normalize();
			double tlight;
			Vector Plight, Nlight;
			int idlight;
			Vector albedo = objects[indexElement]->albedo;
			Ray lightRay(P + 0.001 * N, light);
			double shadow = 1;
			if (sceneIntersection(lightRay, tlight, Plight, Nlight, idlight))
			{

				if (tlight * tlight < distlight)
				{

					shadow = 0;
				}
			}

			color_pixel = shadow * I / ((scene_light - P).norm2() * 4 * M_PI) * albedo / M_PI * std::max(0., dot(light, N));
			Vector randomRayDir = random_cos(N);
			double eps = 10E-2;
			Ray randomRay = Ray(P + eps * N, randomRayDir);

			color_pixel = color_pixel + albedo * getColor(randomRay, bounce - 1);
		}
		return color_pixel;
	}

	std::vector<Object *> objects;
	Vector scene_light;
	double I;
};

class TriangleIndices
{
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group){};
	int vtxi, vtxj, vtxk;
	int uvi, uvj, uvk;
	int ni, nj, nk;
	int group;
};

class BVH
{
public:
	BoundingBox *bbox;
	int i0, i1;
	BVH *left, *right;
};

class TriangleMesh : public Object
{
public:
	~TriangleMesh()
	{
		delete bvh;
	}

	TriangleMesh(const Vector &albedo, bool mirror = false, bool transparent = false) : Object(albedo, mirror, transparent)
	{
		bvh = new BVH;
	};

	void scale_and_translate(double s, const Vector &t)
	{
		for (int i = 0; i < vertices.size(); i++)
		{
			vertices[i] = vertices[i] * s + t;
		}
	}

	bool triangleIntersection(const Ray &ray, const Vector &A, const Vector &B, const Vector &C, double &t, Vector &P, Vector &N)
	{
		Vector e2 = C - A;
		Vector e1 = B - A;
		N = cross(e1, e2);
		double uDotN = dot(ray.direction, N);
		Vector cross1 = cross(A - ray.origin, ray.direction);
		double beta = dot(e2, cross1) / uDotN;
		double gamma = -dot(e1, cross1) / uDotN;
		double alpha = 1 - beta - gamma;
		t = dot(A - ray.origin, N) / uDotN;
		if (beta < 0)
			return false;
		if (gamma < 0)
			return false;
		if (alpha < 0)
			return false;
		if (beta > 1)
			return false;
		if (gamma > 1)
			return false;
		if (alpha > 1)
			return false;
		if (t < 0)
			return false;
		P = ray.origin + t * ray.direction;
		return true;
	}

	bool intersection(const Ray &ray, Vector &P, Vector &N, double &t) override
	{
		double t_inter;
		Vector P_inter, N_inter;
		if (!bvh->bbox->rayBoxIntersection(ray))
			return false;
		bool inter = false;
		std::list<BVH *> BVH_list;

		BVH_list.push_back(bvh);

		while (!BVH_list.empty())
		{
			BVH *currentNode = BVH_list.back();
			BVH_list.pop_back();
			if (currentNode->left && currentNode->bbox->rayBoxIntersection(ray))
			{
				BVH_list.push_back(currentNode->left);
			}
			if (currentNode->right && currentNode->bbox->rayBoxIntersection(ray))
			{
				BVH_list.push_back(currentNode->right);
			}
			if (!currentNode->left)
			{
				for (int i = currentNode->i0; i < currentNode->i1; ++i)
				{
					Vector A = vertices[indices[i].vtxi];
					Vector B = vertices[indices[i].vtxj];
					Vector C = vertices[indices[i].vtxk];
					bool has_inter = triangleIntersection(ray, A, B, C, t_inter, P_inter, N_inter);
					if (has_inter && (t_inter < t))
					{
						inter = true;
						t = t_inter;
						P = P_inter;
						N = N_inter;
					}
				}
			}
		}
		return inter;
	}

	void readOBJ(const char *obj)
	{

		char matfile[255];
		char grp[255];

		FILE *f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f))
		{
			char line[255];
			if (!fgets(line, 255, f))
				break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's')
			{
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ')
			{
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6)
				{
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);
				}
				else
				{
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n')
			{
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't')
			{
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f')
			{
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char *consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9)
				{
					if (i0 < 0)
						t.vtxi = vertices.size() + i0;
					else
						t.vtxi = i0 - 1;
					if (i1 < 0)
						t.vtxj = vertices.size() + i1;
					else
						t.vtxj = i1 - 1;
					if (i2 < 0)
						t.vtxk = vertices.size() + i2;
					else
						t.vtxk = i2 - 1;
					if (j0 < 0)
						t.uvi = uvs.size() + j0;
					else
						t.uvi = j0 - 1;
					if (j1 < 0)
						t.uvj = uvs.size() + j1;
					else
						t.uvj = j1 - 1;
					if (j2 < 0)
						t.uvk = uvs.size() + j2;
					else
						t.uvk = j2 - 1;
					if (k0 < 0)
						t.ni = normals.size() + k0;
					else
						t.ni = k0 - 1;
					if (k1 < 0)
						t.nj = normals.size() + k1;
					else
						t.nj = k1 - 1;
					if (k2 < 0)
						t.nk = normals.size() + k2;
					else
						t.nk = k2 - 1;
					indices.push_back(t);
				}
				else
				{
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6)
					{
						if (i0 < 0)
							t.vtxi = vertices.size() + i0;
						else
							t.vtxi = i0 - 1;
						if (i1 < 0)
							t.vtxj = vertices.size() + i1;
						else
							t.vtxj = i1 - 1;
						if (i2 < 0)
							t.vtxk = vertices.size() + i2;
						else
							t.vtxk = i2 - 1;
						if (j0 < 0)
							t.uvi = uvs.size() + j0;
						else
							t.uvi = j0 - 1;
						if (j1 < 0)
							t.uvj = uvs.size() + j1;
						else
							t.uvj = j1 - 1;
						if (j2 < 0)
							t.uvk = uvs.size() + j2;
						else
							t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else
					{
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3)
						{
							if (i0 < 0)
								t.vtxi = vertices.size() + i0;
							else
								t.vtxi = i0 - 1;
							if (i1 < 0)
								t.vtxj = vertices.size() + i1;
							else
								t.vtxj = i1 - 1;
							if (i2 < 0)
								t.vtxk = vertices.size() + i2;
							else
								t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else
						{
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0)
								t.vtxi = vertices.size() + i0;
							else
								t.vtxi = i0 - 1;
							if (i1 < 0)
								t.vtxj = vertices.size() + i1;
							else
								t.vtxj = i1 - 1;
							if (i2 < 0)
								t.vtxk = vertices.size() + i2;
							else
								t.vtxk = i2 - 1;
							if (k0 < 0)
								t.ni = normals.size() + k0;
							else
								t.ni = k0 - 1;
							if (k1 < 0)
								t.nj = normals.size() + k1;
							else
								t.nj = k1 - 1;
							if (k2 < 0)
								t.nk = normals.size() + k2;
							else
								t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true)
				{
					if (consumedline[0] == '\n')
						break;
					if (consumedline[0] == '\0')
						break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3)
					{
						if (i0 < 0)
							t2.vtxi = vertices.size() + i0;
						else
							t2.vtxi = i0 - 1;
						if (i2 < 0)
							t2.vtxj = vertices.size() + i2;
						else
							t2.vtxj = i2 - 1;
						if (i3 < 0)
							t2.vtxk = vertices.size() + i3;
						else
							t2.vtxk = i3 - 1;
						if (j0 < 0)
							t2.uvi = uvs.size() + j0;
						else
							t2.uvi = j0 - 1;
						if (j2 < 0)
							t2.uvj = uvs.size() + j2;
						else
							t2.uvj = j2 - 1;
						if (j3 < 0)
							t2.uvk = uvs.size() + j3;
						else
							t2.uvk = j3 - 1;
						if (k0 < 0)
							t2.ni = normals.size() + k0;
						else
							t2.ni = k0 - 1;
						if (k2 < 0)
							t2.nj = normals.size() + k2;
						else
							t2.nj = k2 - 1;
						if (k3 < 0)
							t2.nk = normals.size() + k3;
						else
							t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else
					{
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2)
						{
							if (i0 < 0)
								t2.vtxi = vertices.size() + i0;
							else
								t2.vtxi = i0 - 1;
							if (i2 < 0)
								t2.vtxj = vertices.size() + i2;
							else
								t2.vtxj = i2 - 1;
							if (i3 < 0)
								t2.vtxk = vertices.size() + i3;
							else
								t2.vtxk = i3 - 1;
							if (j0 < 0)
								t2.uvi = uvs.size() + j0;
							else
								t2.uvi = j0 - 1;
							if (j2 < 0)
								t2.uvj = uvs.size() + j2;
							else
								t2.uvj = j2 - 1;
							if (j3 < 0)
								t2.uvk = uvs.size() + j3;
							else
								t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else
						{
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2)
							{
								if (i0 < 0)
									t2.vtxi = vertices.size() + i0;
								else
									t2.vtxi = i0 - 1;
								if (i2 < 0)
									t2.vtxj = vertices.size() + i2;
								else
									t2.vtxj = i2 - 1;
								if (i3 < 0)
									t2.vtxk = vertices.size() + i3;
								else
									t2.vtxk = i3 - 1;
								if (k0 < 0)
									t2.ni = normals.size() + k0;
								else
									t2.ni = k0 - 1;
								if (k2 < 0)
									t2.nj = normals.size() + k2;
								else
									t2.nj = k2 - 1;
								if (k3 < 0)
									t2.nk = normals.size() + k3;
								else
									t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else
							{
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1)
								{
									if (i0 < 0)
										t2.vtxi = vertices.size() + i0;
									else
										t2.vtxi = i0 - 1;
									if (i2 < 0)
										t2.vtxj = vertices.size() + i2;
									else
										t2.vtxj = i2 - 1;
									if (i3 < 0)
										t2.vtxk = vertices.size() + i3;
									else
										t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else
								{
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}
			}
		}
		fclose(f);
	}

	BoundingBox *computeBox(int i0, int i1)
	{

		Vector m(1E9, 1E9, 1E9);
		Vector M(-1E9, -1E9, -1E9);
		for (int i = i0; i < i1; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				m[j] = std::min(m[j], std::min(vertices[indices[i].vtxi][j], std::min(vertices[indices[i].vtxj][j], vertices[indices[i].vtxk][j])));
				M[j] = std::max(M[j], std::max(vertices[indices[i].vtxi][j], std::max(vertices[indices[i].vtxj][j], vertices[indices[i].vtxk][j])));
			}
		}
		BoundingBox *res = new BoundingBox;
		res->m = m;
		res->M = M;

		return res;
	}

	void buildBVH(BVH *bvh, int i0, int i1)
	{
		bvh->bbox = computeBox(i0, i1);
		bvh->i0 = i0;
		bvh->i1 = i1;
		bvh->left = NULL;
		bvh->right = NULL;

		Vector diag = bvh->bbox->M - bvh->bbox->m;
		int axis = 0;
		if (diag[1] > diag[0] && diag[1] > diag[2])
		{
			axis = 1;
		}
		else
		{
			if (diag[2] > diag[0] && diag[2] > diag[1])
				axis = 2;
		}

		double splitVal = diag[axis] * 0.5 + bvh->bbox->m[axis];
		int pivot_index = i0;

		std::cout << pivot_index << std::endl;

		for (int i = i0; i < i1; ++i)
		{
			double barycenter = (vertices[indices[i].vtxi][axis] + vertices[indices[i].vtxj][axis] + vertices[indices[i].vtxk][axis]) / 3;
			if (barycenter < splitVal)
			{
				std::swap(indices[i], indices[pivot_index]);
				pivot_index++;
			}
		}

		if ((pivot_index <= i0) || (pivot_index >= i1 - 1) || ((i1 - i0) < 5))
		{
			return;
		}

		std::cout << pivot_index << std::endl;
		std::cout << i1 << std::endl;

		bvh->left = new BVH;
		bvh->right = new BVH;
		buildBVH(bvh->left, i0, pivot_index);
		buildBVH(bvh->right, pivot_index, i1);
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	BVH *bvh;
};

Ray get_light_ray(int i, int j, int H, int W, double fov, Vector Q)
{

	double stdev = 0.5;
	double r1 = uniform(engine);
	double r2 = uniform(engine);
	double x = std::sqrt(-2 * std::log(r1)) * std::cos(2 * M_PI * r2) * stdev;
	double y = std::sqrt(-2 * std::log(r1)) * std::sin(2 * M_PI * r2) * stdev;
	Vector u((j + y) - W / 2 + 0.5, H / 2 - (i + x) - 0.5, -W / (2 * tan(fov / 2)));
	u.normalize();
	return Ray(Q, u);
}

int main()
{

	auto start = clock();

	int W = 512;
	int H = 512;

	double fov = 60 * M_PI / 180;
	Vector Q(0, 0, 55);
	std::vector<unsigned char> image(W * H * 3, 0);
	double I = 3E10;
	int pathsNb = 32;
	Vector light(-10, 20, 40);
	Scene scene(I, light);

	Sphere S1(Vector(0, 0, 0), 10, Vector(1, 1, 1));
	Sphere S2(Vector(-20, 0, 0), 10, Vector(0.7, 0.2, 0.1), true);
	Sphere S3(Vector(20, 0, 0), 10, Vector(0.7, 0.2, 0.1), false, true);

	Sphere Sd(Vector(0, -1000, 0), 990, Vector(0.1, 0.5, 0.1));
	Sphere Su(Vector(0, 1000, 0), 940, Vector(0.2, 0.4, 0.2));
	Sphere Sr(Vector(1000, 0, 0), 940, Vector(0.7, 0.2, 0.5));
	Sphere Sl(Vector(-1000, 0, 0), 940, Vector(0.7, 0.2, 0.5));
	Sphere Sb(Vector(0, 0, -1000), 940, Vector(0.7, 0.2, 0.5));
	Sphere Sf(Vector(0, 0, 1000), 940, Vector(0.9, 0.3, 0.5));

	scene.addToScene(&S1);
	scene.addToScene(&S2);
	scene.addToScene(&S3);

	scene.addToScene(&Sd);
	scene.addToScene(&Su);
	scene.addToScene(&Sr);
	scene.addToScene(&Sl);
	scene.addToScene(&Sb);
	scene.addToScene(&Sf);

	TriangleMesh mesh(Vector(0.8, 0.5, 0.1));
	mesh.readOBJ("cat.obj");
	mesh.scale_and_translate(0.35, Vector(-7, -10, 20));
	mesh.buildBVH(mesh.bvh, 0, mesh.indices.size());
	scene.addToScene(&mesh);

#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < H; ++i)
	{
		for (int j = 0; j < W; ++j)
		{
			Vector color(0, 0, 0);
			for (int k = 0; k < pathsNb; ++k)
			{

				Ray ray = get_light_ray(i, j, H, W, fov, Q);
				color = color + scene.getColor(ray, 5);
			}
			color = color / pathsNb;

			image[(i * W + j) * 3 + 0] = std::min(std::pow(color[0], 0.45), 255.);
			image[(i * W + j) * 3 + 1] = std::min(std::pow(color[1], 0.45), 255.);
			image[(i * W + j) * 3 + 2] = std::min(std::pow(color[2], 0.45), 255.);
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
	std::cout << "finished";
	auto stop = clock();
	double duration_sec = double(stop - start) / CLOCKS_PER_SEC;
	std::cout << " duration: " << duration_sec << " s \n";
	return 0;
}