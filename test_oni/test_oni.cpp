
//
// An example of loading and simulating a rigged model consisting of a collection of rigidbodies and some joints connecting them.
// left mouse and mwheel can select and move individual rigidbodies using on-the-fly positional constraints.
//

#include <map>
#include <strstream>
#include <sstream>
#include <iostream>

#include <mesh.h>
#include <minixml.h>
#include <hull.h>
#include <wingmesh.h>  // just so i can quickly make a box
#include "../testphys/physics.h"
#include <dxwin.h>
#include "ONICamera.h"

template<class T>
std::vector<T*> Addresses(std::vector<T> &a) { std::vector<T*> p; for (auto &e : a) p.push_back(&e); return p; }

//
inline void generatePoints(const uint16_t *depth, const int width, const int height, const float fx, const float fy, const float px, const float py, std::vector<float3> &points) {
	auto halfX = px;
	auto halfY = py;
	auto cX = 1.0f / fx;
	auto cY = 1.0f / fy;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			const auto z = depth[(i * width + j)];
			points[(i * width + j)].x = (j - halfX) * z * cX;
			points[(i * width + j)].y = (i - halfY) * z * cY;
			points[(i * width + j)].z = z;
		}
	}
}

//must be a power of two
template <int numVoxels>
void voxelSubsample(const std::vector<float3> & input, std::vector<float3> &output, float voxelSize, int minVoxelNum)
{
	struct VoxelCandidate {
		int3 pos = { 0, 0, 0 };
		float3 sum = { 0, 0, 0 };
		int cnt = 0;
	};	typedef unsigned int uint;
	static VoxelCandidate cand[numVoxels];
	memset(&cand, 0, numVoxels*sizeof(VoxelCandidate));
	auto voxelMask = numVoxels - 1;
	auto iVS = 1.0f / voxelSize; //inverse voxel size
	//auto hashCoeff = uint3(3863,2029,49139); //good for mm
	//auto hashCoeff = uint3(54851, 11909, 24781); //good for m
	int3 hashCoeff = { 56599, 11399, 33359 };
	auto collisions = 0;
	for (const auto & pt : input) {
		int3 intPos(static_cast<uint>(std::floor(pt.x*iVS)), static_cast<uint>(std::floor(pt.y*iVS)), static_cast<uint>(std::floor(pt.z*iVS)));
		unsigned int hash = dot(hashCoeff, intPos);
		hash &= voxelMask;
		auto & bucket = cand[hash];
		if (bucket.cnt && bucket.pos != intPos) { //flush on collision
			output.emplace_back(bucket.sum / static_cast<float>(bucket.cnt));
			bucket.cnt = 0;
			bucket.sum = { 0, 0, 0 };
			collisions++;
		}
		if (!bucket.cnt) {
			bucket.pos = intPos;
		}
		bucket.sum += pt;
		bucket.cnt++;
	}
	for (const auto &pt : cand)
	{
		if (pt.cnt >= minVoxelNum)
			output.emplace_back(pt.sum / static_cast<float>(pt.cnt));
	}
	//std::cout << hashCoeff[1] << "\t" << collisions << std::endl;
}
//must be a power of two
template <int numVoxels>
std::vector<float3> voxelSubsample(const std::vector<float3> & input, float voxelSize, int minVoxelNum)
{
	std::vector<float3> output;
	voxelSubsample<numVoxels>(input, output, voxelSize, minVoxelNum);
	return output;
}

float4 closestplane(const std::vector<float4> &planes, const float3 &v, const float3 &n)
{
	float4 r = { 0, 0, 0, -std::numeric_limits<float>::max() };
	for (auto p : planes)
	{
		if (dot(n, p.xyz()) < 0)  // to filter non-camera facing planes
			continue;
		if (dot(float4(v, 1), p) > dot(float4(v, 1), r))
			r = p;
	}
	return r;
}


struct Joint
{
	int rbi0;
	int rbi1;
	float3 p0;
	float3 p1;
	float3 jointlimitmin;
	float3 jointlimitmax;
};

std::vector<float4> Planes(const std::vector<float3> &verts, const std::vector<int3> &tris) { std::vector<float4> planes; for (auto &t : tris) planes.push_back(PolyPlane({ verts[t[0]], verts[t[1]], verts[t[2]] }));  return planes; }

int main(int argc, const char *argv[]) try
{
		ONICamera cam;
		auto quit = !cam.Start();
		auto fx = cam.getFx();
		auto fy = cam.getFy();
		auto px = cam.getPx();
		auto py = cam.getPy();
		auto width = cam.getXDim();
		auto height = cam.getYDim();
		std::vector<float3> points(width * height);
	
	std::vector<Joint> joints;
	std::vector<RigidBody> rbs;
	std::map<std::string, unsigned int> rbindex;
	auto xml = XMLParseFile("./default_hand.chr");  // replace string with whatever model you want to test.   uses xml variation of John Ratcliff's easy mesh (ezm) file format.
	auto const &skx = xml.child("model").child("skeleton");
	for (auto const &b : skx.children)
	{
		rbindex[b.attribute("name")] = rbs.size();
		auto verts = ArrayImport<float3>(b.child("verts").body);
		auto tris = calchull(verts, verts.size());
		float3 pos = StringTo<float3>(b.attribute("position"));
		int parent = (b.hasAttribute("parent")) ? (int)rbindex[b.attribute("parent")] : -1;
		rbs.push_back(RigidBody({ Shape(verts, tris) }, pos + ((parent >= 0) ? rbs[parent].position - rbs[parent].com : float3(0, 0, 0))));
		if (parent >= 0)
			joints.push_back({ parent, (int)rbs.size() - 1, pos, float3(0, 0, 0), StringTo<float3>(b.child("jointlimitmin").body), StringTo<float3>(b.child("jointlimitmax").body) });
	}
	rbscalemass(&rbs[0], 3.0f);
	rbscalemass(&rbs[1], 5.0f);

	DXWin mywin("DX testing articulated rigged model");
	//OVRWin mywin("VR testing articulated rigged model");
	std::vector<Mesh> meshes;
	for (auto &rb : rbs)
	{
		meshes.push_back(MeshSmoothish(rb.shapes[0].verts, rb.shapes[0].tris)); //  1 shape each is known
		rb.damping = 0.8f;
		//rb.gravscale = 0;
	}
	for (auto &joint : joints)
	{
		rbs[joint.rbi0].ignore.push_back(&rbs[joint.rbi1]);
		rbs[joint.rbi1].ignore.push_back(&rbs[joint.rbi0]);
		joint.p0 -= rbs[joint.rbi0].com;
		joint.p1 -= rbs[joint.rbi1].com;
	}
	for (auto &ja : joints) for (auto &jb : joints) if (ja.rbi0 == jb.rbi0 && ja.rbi1 != jb.rbi1)  // ignore siblings 
	{
		rbs[ja.rbi1].ignore.push_back(&rbs[jb.rbi1]);
		rbs[jb.rbi1].ignore.push_back(&rbs[ja.rbi1]);
	}
	for (auto &ja : joints) for (auto &jb : joints) if (ja.rbi1 == jb.rbi0)  // ignore grandparents 
	{
		rbs[ja.rbi0].ignore.push_back(&rbs[jb.rbi1]);
		rbs[jb.rbi1].ignore.push_back(&rbs[ja.rbi0]);
	}

	std::vector<float3> groundpoints = { { -5.0f, -5.0f, -5.0f }, { 5.0f, -5.0f, -5.0f }, { 5.0f, 10.0f, -5.0f }, { -5.0f, 10.0f, -5.0f }, { -5.0f, -5.0f, -10.0f }, { 5.0f, -5.0f, -10.0f }, { 5.0f, 10.0f, -10.0f }, { -5.0f, 10.0f, -10.0f } };
	Mesh ground = MeshSmoothish(groundpoints, { { 0, 1, 2 }, { 2, 3, 0 } });
	ground.hack = { 1, 1, 0, 1 };
	WingMesh cube_wm = WingMeshCube(0.025f);
	auto mesh_cube = MeshFlatShadeTex(cube_wm.verts, WingMeshTris(cube_wm));
	mesh_cube.hack = { 0, 1, 0, 1 };

	Pose camera = { { 0, -10, 0 }, normalize(float4(1, 0, 0, 1)) };
	RigidBody *selected = NULL;
	float3 spoint = camera * float3(0, 0, -10);
	float3 rbpoint;

	struct Pin{ float3 w; RigidBody* rb; float3 p; };
	std::vector<Pin> pins;

	mywin.keyboardfunc = [&](int key, int, int)
	{
		if (key == 'g') for (auto &rb : rbs) rb.gravscale = 1.0f - rb.gravscale;
		if (key == 'p' && selected)
			Append<Pin>(pins, { spoint, selected, rbpoint });
	};

	while (mywin.WindowUp() && !quit)
	{
		quit = !cam.syncNext();
		memset(points.data(), 0, sizeof(float3) * width * height);
		generatePoints(cam.getDepth(), width, height, fx, fy, px, py, points);
		std::vector<float3> segmentedPoints;
		auto minE = 32768;
		auto maxE = 0;
		for (const auto pt : points){
			if (pt.z > 0){
				minE = std::min<int>(minE, (int)pt.z);
			}
		}
		maxE = minE + 200;
		for (const auto pt : points){
			if (pt.z > 0.01 && pt.z < maxE){
				segmentedPoints.push_back((1 / 100.0f)*pt);
			}
		}
		auto depthdata = voxelSubsample<1024>(segmentedPoints, 0.25f, 1);
		std::vector<float3> originalData; // generated pointcloud 
		float3 ray = qrot(camera.orientation, normalize(mywin.MouseVector));
		if (!selected)
		{
			for (auto &rb : rbs)
			{
				float3 v1 = camera.position + ray*100.0f;
				if (auto h = ConvexHitCheck(Planes(rb.shapes[0].verts, rb.shapes[0].tris), rb.pose(), camera.position, v1))
				{
					v1 = h.impact;
					selected = &rb;
					spoint = h.impact;
					rbpoint = rb.pose().Inverse()*h.impact;
				}
			}
		}
		spoint = camera.position + ray * magnitude(spoint - camera.position)*powf(1.025f, (float)mywin.mousewheel);
		mesh_cube.pose.position = spoint;
		if (!mywin.MouseState)
			selected = NULL;

		std::vector<LimitAngular> angulars;
		std::vector<LimitLinear>  linears;
		for (auto const &joint : joints)
		{
			Append(linears, ConstrainPositionNailed(&rbs[joint.rbi0], joint.p0, &rbs[joint.rbi1], joint.p1));
			Append(angulars, ConstrainAngularRange(&rbs[joint.rbi0], &rbs[joint.rbi1], { 0, 0, 0, 1 }, joint.jointlimitmin, joint.jointlimitmax));
		}
		if (selected)
			Append(linears, ConstrainPositionNailed(NULL, spoint, selected, rbpoint));
		for (auto &p : pins)
			Append(linears, ConstrainPositionNailed(NULL, p.w, p.rb, p.p));
		std::vector<std::pair<float3, float3>> match;
		std::vector<std::pair<RigidBody*, std::vector<float4>>> cached;
	
		for (auto &rb : rbs)
		{
			auto rb_planes = Planes(rb.shapes[0].verts, rb.shapes[0].tris);
			for (auto & rbp : rb_planes) {
				rbp = rb.pose().TransformPlane(rbp);
			}
			cached.push_back({ &rb, rb_planes });
		}
		for (auto &pt : depthdata) {
			pt = float3(-pt.x, pt.z -5.0, -pt.y + 1.0f);
		}
		for (auto &p0 : depthdata) {
			float closest = INT_MAX;
			std::pair<RigidBody*, std::vector<float4>> *rbcB = nullptr;
			for (auto &rbc : cached)
			{
				auto rb = rbc.first;
				for (const auto vert : rb->shapes[0].verts){
					float dist = mag2(rb->pose()*vert - p0);
					if (dist < closest) {
						closest = dist;
						rbcB = &rbc;
					}
				}
			}

			auto planes = rbcB->second;
			auto plane = closestplane(planes, p0, { 0, 0, 0 });  // could pass normal direction of -p0 to avoid backside planes

			//if (plane.w < 0)
			//	plane.w = std::min(0.0f, plane.w + dot(plane.xyz(), normalize(p0))*0.2f);  // small hack here (may add jitter)!! add thickness if we are using a backside plane
			auto p1w = p0 - plane.xyz()*dot(plane, float4(p0, 1));               // p1 is on the plane
			match.push_back(std::pair<float3, float3>(p0, p1w));

			linears.push_back(ConstrainAlongDirection(NULL, p0, rbcB->first, rbcB->first->pose().Inverse()*p1w, plane.xyz(), -50, 50));
		}
		PhysicsUpdate(Addresses(rbs), linears, angulars, { &groundpoints });

		for (unsigned int i = 0; i < rbs.size(); i++)
		{
			meshes[i].pose = rbs[i].pose();
		}
		std::vector<Mesh> dpt;
		std::vector<Mesh*> mshptr = { &ground, &mesh_cube };
		for (const auto &pt : depthdata) {
			Mesh newCube = mesh_cube;
			newCube.pose.position = pt;
			dpt.push_back(newCube);
		}
		for (auto &pt : dpt) {
			mshptr.push_back(&pt);
		}
		mywin.RenderScene(camera, Append(Addresses(meshes), mshptr));
	}
}
catch (std::exception e)
{
	MessageBoxA(GetActiveWindow(), e.what(), "FAIL", 0);
	return -1;
}