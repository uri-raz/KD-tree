#include "stdafx.h"

#include <ctime>
#include <random>
#include <iostream>

#include "kd tree.h"

void Cube(std::vector<kd_tree<3>::kd_point> &Points, const kd_tree<3>::kd_point &center, float size)
{
	Points.push_back(center);

	if (size > 1.0f)

		for (unsigned i = 0; i < 8; i++)
		{
			float mul[3] = { i & 1 ? -1.0f : 1.0f ,
							 i & 2 ? -1.0f : 1.0f ,
							 i & 4 ? -1.0f : 1.0f };

			kd_tree<3>::kd_point corner = { center[0] + mul[0] * size,
										    center[1] + mul[1] * size,
										    center[2] + mul[2] * size };

			Cube(Points, corner, size / 4.0f);
		}
}

void performance(const kd_tree<3> &Tree)
{
	size_t errors = 0;

	// Prep sample points for testing

	kd_box<3> minmax = *Tree.boundingBox();

	std::default_random_engine generator;

	std::uniform_real_distribution<float> distributionX(minmax.first[0] - 1.5f, minmax.second[0] + 1.5f);
	std::uniform_real_distribution<float> distributionY(minmax.first[1] - 1.5f, minmax.second[1] + 1.5f);
	std::uniform_real_distribution<float> distributionZ(minmax.first[2] - 1.5f, minmax.second[2] + 1.5f);

	const unsigned sampleSize = 10000;
	std::vector<kd_tree<3>::kd_point> samplePoints;
	samplePoints.resize(sampleSize);

	for (unsigned index = 0; index < sampleSize; index++)
		samplePoints[index] = kd_tree<3>::kd_point({ distributionX(generator),
													 distributionY(generator),
													 distributionZ(generator) });

	// Test KNearestNeighbors

	const unsigned setSize = 16;

	auto startTime = clock();

	for (const auto& Point: samplePoints)
	{
		auto nearPoints = Tree.KNearestNeighbors(Point, setSize);

		if (nearPoints.has_value() && nearPoints.value().size() != setSize)
			errors++;
	}

	auto endTime = clock();
	auto runTime = (endTime - startTime) / static_cast<double>(CLOCKS_PER_SEC);

	std::cout << "Did " << sampleSize << " K nearest neighbors searches in " << runTime << " seconds.\n";

	if (errors != 0)
		std::cout << "Had " << errors << " search errors.\n";

	errors = 0;

	startTime = clock();

	for (const auto& Point : samplePoints)
		if (!Tree.nearestNeighbor(Point))
			errors++;

	endTime = clock();
	runTime = (endTime - startTime) / static_cast<double>(CLOCKS_PER_SEC);

	if (errors != 0)
		std::cout << "Had " << errors << " search errors.\n";

	std::cout << "Did " << sampleSize << " nearest neighbor searches in " << runTime << " seconds.\n";

	// Test K Nearest Neighbors for 1 point

	errors = 0;
	startTime = clock();


	for (const auto& Point : samplePoints)
		if (Tree.KNearestNeighbors(Point, 1).value().size() != 1)
			errors++;

	endTime = clock();
	runTime = (endTime - startTime) / static_cast<double>(CLOCKS_PER_SEC);

	if (errors != 0)
		std::cout << "Had " << errors << " search errors.\n";

	std::cout << "Did " << sampleSize << " single nearest neighbor searches in " << runTime << " seconds.\n";

	std::cout << std::endl;
}

int _tmain(int argc, _TCHAR* argv[])
{
	// Code below is toy test

	kd_tree<3> Tree;
	kd_tree<3>::kd_point Point = { 0.5f, 0.5f, 0.5f },
						 MinCorner = { 0.0f, 0.0f, 0.0f },
						 MaxCorner = { 2.0f, 2.0f, 2.0f };

	std::vector<kd_tree<3>::kd_point> Points;
	std::optional<std::vector<kd_tree<3>::kd_point>> nearPoints;

	// Create a tree

	Points = { kd_tree<3>::kd_point({ 0.0f, 0.0f, 0.0f }),
			   kd_tree<3>::kd_point({ 0.0f, 0.0f, 1.0f }),
			   kd_tree<3>::kd_point({ 0.0f, 1.0f, 0.0f }),
			   kd_tree<3>::kd_point({ 0.0f, 1.0f, 1.0f }),
			   kd_tree<3>::kd_point({ 1.0f, 0.0f, 0.0f }),
			   kd_tree<3>::kd_point({ 1.0f, 0.0f, 1.0f }),
			   kd_tree<3>::kd_point({ 1.0f, 1.0f, 0.0f }),
			   kd_tree<3>::kd_point({ 1.0f, 1.0f, 1.0f }),
			   kd_tree<3>::kd_point({ 0.0f, 0.0f, 2.0f }),
			   kd_tree<3>::kd_point({ 0.0f, 2.0f, 0.0f }),
			   kd_tree<3>::kd_point({ 0.0f, 2.0f, 2.0f }),
			   kd_tree<3>::kd_point({ 2.0f, 0.0f, 0.0f }),
			   kd_tree<3>::kd_point({ 2.0f, 0.0f, 2.0f }),
			   kd_tree<3>::kd_point({ 2.0f, 2.0f, 0.0f }),
			   kd_tree<3>::kd_point({ 2.0f, 2.0f, 2.0f }) };

	for (auto const &p : Tree)
		std::cout << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
	std::cout.flush();
	
	Tree.insert(Points);
	Tree.PrintTree();

	// Call tree height & nodes count methods

	std::cout << "Tree height is " << Tree.TreeHeight() << '\n';
	std::cout << "There are " << Tree.nodeCount() << " points in the tree.\n";

	// Use an iterator to print all points in tree

	std::cout << std::endl;
	std::cout << "Points in cloud are:\n";

	for (const auto& point : Tree)
		std::cout << point[0] << " " << point[1] << " " << point[2] << '\n';

	std::cout << std::endl;

	auto nearPoint = Tree.nearestNeighbor(Point);

	if (nearPoint)
		std::cout << "Closest point to ("
		<< Point[0] << "," << Point[1] << "," << Point[2] << ") is ("
		<< nearPoint.value()[0] << "," << nearPoint.value()[1] << "," << nearPoint.value()[2] << ")" << '\n';
	else
		std::cout << "Error: no nearest point found!\n";

	nearPoints = Tree.KNearestNeighbors(Point, 8);

	std::cout << std::endl;
	std::cout << "8 closest points to ("
		<< Point[0] << "," << Point[1] << "," << Point[2] << ") are\n";

	{
		auto const& nPoints = nearPoints.value();
		for (size_t i = 0; i < nPoints.size(); i++)
			std::cout << " (" << nPoints[i][0] << "," << nPoints[i][1] << "," << nPoints[i][2] << ")" << std::endl;
	}

	Tree.pointsInBox(kd_box<3>(MinCorner, MaxCorner), Points);

	std::cout << std::endl;
	std::cout << "Points in bounding box ";
	std::cout << " (" << MinCorner[0] << "," << MinCorner[1] << "," << MinCorner[2] << ") ->";
	std::cout << " (" << MaxCorner[0] << "," << MaxCorner[1] << "," << MaxCorner[2] << ") are\n";

	for (size_t i = 0; i < Points.size(); i++)
		std::cout << " (" << Points[i][0] << "," << Points[i][1] << "," << Points[i][2] << ")\n";

	std::cout << std::endl;
	std::cout << "Tree dump is" << std::endl;
	Tree.PrintTree();

	std::cout << std::endl;

	Tree.clear();
	Points.clear();

	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 0.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 0.0f, 1.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 1.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 1.0f, 1.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 1.0f, 0.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 1.0f, 0.0f, 1.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 1.0f, 1.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 1.0f, 1.0f, 1.0f }));
	Tree.insert(Points);
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 0.0f, 0.0f, 0.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 0.0f, 0.0f, 1.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 1.0f, 0.0f, 0.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 1.0f, 0.0f, 1.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 0.0f, 1.0f, 1.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 1.0f, 1.0f, 1.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 0.0f, 1.0f, 0.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 1.0f, 1.0f, 0.0f }));
	Tree.PrintTree(); std::cout << std::endl;

	Tree.clear();

	Tree.insert(kd_tree<3>::kd_point({ 1.0f, 1.0f, 0.0f }));
	Tree.insert(kd_tree<3>::kd_point({ 1.0f, 1.0f, 1.0f }));
	Tree.insert(kd_tree<3>::kd_point({ 1.0f, 0.0f, 0.0f }));
	Tree.insert(kd_tree<3>::kd_point({ 1.0f, 0.0f, 1.0f }));
	Tree.insert(kd_tree<3>::kd_point({ 0.0f, 1.0f, 0.0f }));
	Tree.insert(kd_tree<3>::kd_point({ 0.0f, 1.0f, 1.0f }));
	Tree.insert(kd_tree<3>::kd_point({ 0.0f, 0.0f, 0.0f }));
	Tree.insert(kd_tree<3>::kd_point({ 0.0f, 0.0f, 1.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 0.0f, 0.0f, 0.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 0.0f, 0.0f, 1.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 1.0f, 0.0f, 0.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 1.0f, 0.0f, 1.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 0.0f, 1.0f, 1.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 1.0f, 1.0f, 1.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 0.0f, 1.0f, 0.0f }));
	Tree.PrintTree(); std::cout << std::endl;
	Tree.erase(kd_tree<3>::kd_point({ 1.0f, 1.0f, 0.0f }));
	Tree.PrintTree(); std::cout << std::endl;

	Tree.clear();
	Points.clear();

	// -----------------------

	// Code below is unit test

	time_t startTime = clock();

	std::vector<kd_tree<3>::kd_point> PointsB;


	if (Tree.TreeHeight() == 0)
		std::cout << "OK: height of empty tree is zero." << std::endl;
	else
		std::cout << "Error: height of empty tree is not zero." << std::endl;

	Tree.pointsInBox(kd_box<3>(kd_tree<3>::kd_point({ -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity() }),
							   kd_tree<3>::kd_point({ std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity() })),
							   Points);

	if (Points.size() == 0)
		std::cout << "OK: searched for all points in an empty tree, and found none." << std::endl;
	else
		std::cout << "Error: searched for all points in an empty tree, and found " << Points.size() << std::endl;

	Tree.insert(kd_tree<3>::kd_point({ 0, 0, 0 }));

	Tree.pointsInBox(kd_box<3>(kd_tree<3>::kd_point({ -1.0, -1.0, -1.0 }), kd_tree<3>::kd_point({ 1.0, 1.0, 1.0 })), Points);

	if (Points.size() == 1)
		std::cout << "OK: searched for points in a one points tree, and found it." << std::endl;
	else
		std::cout << "Error: searched for points in a one points tree, and found " << Points.size() << std::endl;

	Tree.pointsInBox(kd_box<3>(kd_tree<3>::kd_point({ -1.0, -1.0, -1.0 }), kd_tree<3>::kd_point({ -.5, -.5, -.5 })), Points);

	if (Points.size() == 0)
		std::cout << "OK: searched for points using non-crossing box, and found nothing." << std::endl;
	else
		std::cout << "Error: searched for points using non-crossing box, and found " << Points.size() << " points." << std::endl;

	Tree.clear();
	Points.clear();


	Tree.insert(Points);

	if (Tree.TreeHeight() == 0)
		std::cout << "OK: height of tree made from no points is zero." << std::endl;
	else
		std::cout << "Error: height of tree made from no points is non-zero." << std::endl;

	
	Point = { 1.0f, 2.0f, 3.0f };
	Points.push_back(Point);
	Tree.insert(Points);

	if (Tree.TreeHeight() == 1)
		std::cout << "OK: height of tree with one point is one." << std::endl;
	else
		std::cout << "Error: height of tree with one point is not one." << std::endl;

	
	auto Iter = Tree.begin();

	if (*Iter == Point)
		std::cout << "OK: basic iterator test passed." << std::endl;
	else
		std::cout << "Error: basic iterator test failed." << std::endl;
	

	Points.clear();
	Tree.clear();
	
	Points.push_back(Point);
	Points.push_back(kd_tree<3>::kd_point({ std::numeric_limits<float>::quiet_NaN(), 3.5f, -2.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 3.5f, std::numeric_limits<float>::infinity(), -2.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 3.5f, -2.0f, -std::numeric_limits<float>::infinity() }));

	Tree.insert(Points);

	if (Tree.TreeHeight() == 1)
		std::cout << "OK: invalid points filtered." << std::endl;
	else
		std::cout << "Error: invalid points improperly filtered." << std::endl;

	nearPoints = Tree.KNearestNeighbors(kd_tree<3>::kd_point({ 0.0f, 0.0f, 0.0f }), 27);
	if (nearPoints.has_value() && (*nearPoints).size() == 1)
		std::cout << "OK: found 1 nearest neighbor, as expected." << std::endl;
	else
		std::cout << "Error: found " << (*nearPoints).size() << " nearest points in a tree with a single point." << std::endl;

	Points.clear();
	Tree.clear();

	for (auto x = -2.0f; x < 3.0f; x++)
	for (auto y = -2.0f; y < 3.0f; y++)
	for (auto z = -2.0f; z < 3.0f; z++)
		Points.push_back(kd_tree<3>::kd_point({ x, y, z }));

	Tree.insert(Points);

	const kd_tree<3>::kd_point searchPoint = { 0.0f, 0.0f, 0.0f };

	nearPoints = Tree.KNearestNeighbors(searchPoint, 27);
	if (nearPoints.has_value() && (*nearPoints).size() == 27)
	{
		std::sort(Points.begin(), Points.end(),
			      [searchPoint](const kd_tree<3>::kd_point &A, const kd_tree<3>::kd_point &B)
					{return kd_tree<3>::DistanceSq(searchPoint, A) < kd_tree<3>::DistanceSq(searchPoint, B); });
															// Get the 27 closest points (using Distance) at the beginning
		std::sort(Points.begin(), Points.begin() + 27);		// Order those 27 points 'naturally' (using kd_point operator<)

		auto& nPoints = nearPoints.value();

		std::sort(nPoints.begin(), nPoints.end());			// Order the nearest points 'naturally' (using kd_point operator<)

		bool equal = true;

		for (unsigned i = 0; i < 27; i++)
			if (Points[i] != nPoints[i])
			{
				equal = false;
				break;
			}

		if (equal)
			std::cout << "OK: Correct k nearest points found" << std::endl;
		else
			std::cout << "Error: incorrect k nearest points found" << std::endl;
	}
	else
		std::cout << "Error: found incorrect number of nearest points in a non-empty tree." << std::endl;


	for (const auto& point : Tree)
		PointsB.push_back(point);

	if (PointsB.size() == Points.size())
	{
		sort(Points.begin(), Points.end());
		sort(PointsB.begin(), PointsB.end());

		bool equal = true;
		size_t i;
		
		for (i = 0; i < Points.size(); i++)
			if (Points[i] != PointsB[i])
			{
				equal = false;
				break;
			}

			if (equal)
				std::cout << "OK: iterator got all tree points correctly." << std::endl;
			else
				std::cout << "Error: iterator got an incorrect point from tree, see index " << i << std::endl;
	}
	else
		std::cout << "Error: inserted " << Points.size() << " points to tree, and iterator got " << PointsB.size() << std::endl;

	Tree.clear();
	Points.clear();
	PointsB.clear();


	Cube(Points, kd_tree<3>::kd_point({ 0.0f, 0.0f, 0.0f }), 16.0f);
	Tree.insert(Points);

	auto nodeCount = Tree.nodeCount();

	if (nodeCount != Points.size())
		std::cout << "Error: cubical tree should have " << Points.size() << " nodes, and has " << nodeCount << std::endl;
	else
		std::cout << "OK: cubical tree has the expected number of nodes." << std::endl;

	kd_box<3> refBox{ kd_tree<3>::kd_point{-20.0f, -20.0f, -20.0f}, kd_tree<3>::kd_point{ 20.0f, 20.0f, 20.0f } };

	if (Tree.boundingBox() == refBox)
		std::cout << "OK: bounding box is correct." << std::endl;
	else
		std::cout << "Error: bounding box is incorrect." << std::endl;

	{
		std::default_random_engine generator;
		std::uniform_real_distribution<float> distribution(-10000, 10000);

		for (unsigned index = 0; index < 1024; index++)
		{
			const auto randPoint = kd_tree<3>::kd_point({ distribution(generator),
														  distribution(generator),
														  distribution(generator) });

			if (auto point = Tree.nearestNeighbor(randPoint))
			{
				auto nearPoint = Points[0];
				auto minDistance = kd_tree<3>::Distance(randPoint, Points[0]);

				for (const auto& iPoint : Points)
					if (kd_tree<3>::Distance(randPoint, iPoint) < minDistance)
					{
						nearPoint = iPoint;
						minDistance = kd_tree<3>::Distance(randPoint, nearPoint);
					}

				if (kd_tree<3>::Distance(randPoint, nearPoint) != kd_tree<3>::Distance(randPoint, *point))
				{
					std::cout << "Error: did not find nearest point for (" << randPoint[0] << ","
																		   << randPoint[1] << ","
																		   << randPoint[2] << ")" << std::endl;
					std::cout << "Found (" << (*point)[0] << "," << (*point)[1] << "," << (*point)[2] << ")" <<
						" " << kd_tree<3>::Distance(randPoint, *point) << std::endl;
					std::cout << "Actual (" << nearPoint[0] << "," << nearPoint[1] << "," << nearPoint[2] << ")" <<
						" " << kd_tree<3>::Distance(randPoint, nearPoint) << std::endl;
				}
			}
			else
				std::cout << "Error: did not find nearest point for (" << randPoint[0] << ","
																	   << randPoint[1] << ","
																	   << randPoint[2] << ")" << std::endl;
		}
	}

	Points.clear();
	Tree.clear();

	if (Tree.nearestNeighbor(kd_tree<3>::kd_point({ 0, 0, 0 })))
		std::cout << "Error: searched for nearest point in an empty tree, and found a point!" << std::endl;
	else
		std::cout << "OK: searched for nearest point in an empty tree, and found nothing." << std::endl;

	Tree.pointsInBox(kd_box<3>(kd_tree<3>::kd_point({ -1.0, -1.0, -1.0 }), kd_tree<3>::kd_point({ 1.0, 1.0, 1.0 })), Points);

	if (Points.size() != 0)
		std::cout << "Error: searched for points in an empty tree, and found " << Points.size() << " points!" << std::endl;
	else
		std::cout << "OK: searched for points in an empty tree, and found none." << std::endl;

	Tree.clear();

	for (float i = -2.0; i < 3.0; i++)
	for (float j = -2.0; j < 3.0; j++)
	for (float k = -2.0; k < 3.0; k++)
		Points.push_back(kd_tree<3>::kd_point({ i, j, k }));

	Tree.insert(Points);

	if (Tree.nodeCount() == 125)
		std::cout << "OK: inserted 125 points, and found 125 nodes in tree." << std::endl;
	else
		std::cout << "Error: inserted 125 points, and found " << Tree.nodeCount() << " nodes in tree." << std::endl;

	Tree.pointsInBox(kd_box<3>(kd_box<3>::kd_point({ -2.0, -2.0, -2.0 }), kd_box<3>::kd_point({ 2.0, 2.0, 2.0 })), Points);

	if (Points.size() == pow(5, 3))
		std::cout << "OK: inserted 125 points, and found them all." << std::endl;
	else
		std::cout << "Error: searched for 125 points, and found " << Points.size() << std::endl;

	Tree.pointsInBox(kd_box<3>(kd_box<3>::kd_point({ 0.1f, 0.1f, 0.1f }), kd_box<3>::kd_point({ 2.1f, 2.1f, 2.1f })), Points);

	if (Points.size() == pow(2, 3))
		std::cout << "OK: quadrant search yielded 8 points, as expected." << std::endl;
	else
		std::cout << "Error: searched for 8 points, and found " << Points.size() << std::endl;

	Tree.clear();

	if (Tree.nodeCount(true) == 0)
		std::cout << "OK: cleared tree has no nodes." << std::endl;
	else
		std::cout << "Error: cleared tree has " << Tree.nodeCount(true) << " nodes." << std::endl;

	Points.clear();
	for (unsigned h = 0; h < 2; h++)
	for (float i = -9.0; i < 10.0; i++)
	for (float j = -9.0; j < 10.0; j++)
	for (float k = -9.0; k < 10.0; k++)
		Points.push_back(kd_tree<3>::kd_point({ i, j, k }));

	Tree.insert(Points);
	Tree.pointsInBox(kd_box<3>(kd_tree<3>::kd_point({ -2.0, -2.0, -2.0 }), kd_tree<3>::kd_point({ 2.0, 2.0, 2.0 })), Points);

	if (Points.size() == pow(5, 3))
		std::cout << "OK: inserted 125 points in triplicity, and only one copy of each was inserted." << std::endl;
	else
		std::cout << "Error: inserted 125 points in triplicity, and found " << Points.size() << std::endl;

	Tree.clear();
	Points.clear();


	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution(-10, 10);

	for (unsigned index = 0; index < 1000; index++)
		Points.push_back(kd_tree<3>::kd_point({ distribution(generator), distribution(generator), distribution(generator) }));

	Tree.insert(Points);

	bool insertedPointsFlag = true;
	for (size_t i = 0; i < Points.size(); i++)
	{
		if (auto point = Tree.nearestNeighbor(Points[i]))
		{
			if (Points[i] != *point)
			{
				insertedPointsFlag = false;

				std::cout << "Error: nearest neighbor of inserted point ("
					<< Points[i][0] << "," << Points[i][1] << "," << Points[i][2] << ") is ("
					<< (*point)[0] << "," << (*point)[1] << "," << (*point)[2] << ")" << std::endl;
			}
		}
		else
		{
			insertedPointsFlag = false;

			std::cout << "Error: failed to find inserted point ("
				<< Points[i][0] << "," << Points[i][1] << "," << Points[i][2] << ")" << std::endl;
		}
	}

	if (insertedPointsFlag)
		std::cout << "OK: nearest neighbor of inserted point is itself." << std::endl;

	bool nearestSearchFlag = true;

	for (unsigned index = 0; index < 1000; index++)
	{
		const auto randPoint = kd_tree<3>::kd_point({ distribution(generator),
													  distribution(generator),
													  distribution(generator) });

		if (const auto point = Tree.nearestNeighbor(randPoint))
		{
			nearPoint = Points[0];
			auto minDistance = kd_tree<3>::Distance(randPoint, Points[0]);

			for (const auto& iPoint : Points)
				if (kd_tree<3>::Distance(randPoint, iPoint) < minDistance)
				{
					nearPoint = iPoint;
					minDistance = kd_tree<3>::Distance(randPoint, iPoint);
				}

			if (kd_tree<3>::Distance(randPoint, *point) != kd_tree<3>::Distance(randPoint, *nearPoint))
			{
				std::cout << "Error: did not find nearest point for (" << randPoint[0] << ","
																	   << randPoint[1] << ","
																	   << randPoint[2] << ")" << std::endl;
				std::cout << "Found (" << (*point)[0] << "," << (*point)[1] << "," << (*point)[2] << ")" <<
					" " << kd_tree<3>::Distance(randPoint, *point) << std::endl;
				std::cout << "Actual (" << (*nearPoint)[0] << "," << (*nearPoint)[1] << "," << (*nearPoint)[2] << ")" <<
					" " << kd_tree<3>::Distance(randPoint, *nearPoint) << std::endl;
			}
		}
		else
		{
			nearestSearchFlag = false;

			std::cout << "Error: did not find nearest point for (" << randPoint[0] << ","
																   << randPoint[1] << ","
																   << randPoint[2] << ")" << std::endl;
		}
	}

	if (nearestSearchFlag)
		std::cout << "OK: nearest neighbor search is in order." << std::endl;

	std::cout << std::endl;

	Tree.clear();
	Points.clear();

	time_t endTime = clock();
	auto runTime = (endTime - startTime) / static_cast<double>(CLOCKS_PER_SEC);
	std::cout << "Run time is " << runTime << " seconds." << std::endl;

	std::cout << std::endl;

	// ------------------------------

	// Code below is performance test

	std::cout << "Reading file" << std::endl;

	FILE *tooth = fopen("D:\\Uri Raz\\C++ Projects\\Tooth filtered from 116 photos.txt", "r");

	auto status = fscanf(tooth, "%f %f %f", &Point[0], &Point[1], &Point[2]);

	while (status != EOF)
	{
		Points.push_back(Point);
		status = fscanf(tooth, "%f %f %f", &Point[0], &Point[1], &Point[2]);
	}

	fclose(tooth);

	std::cout << std::endl;

	// Test tree built the inefficient way
	/*
	std::cout << "Phase I" << std::endl;
	std::cout << "Building tree" << std::endl;

	startTime = clock();
	for (size_t i = 0; i < Points.size(); i++)
		Tree.insert(Points[i]);
	endTime = clock();

	runTime = (endTime - startTime) / (float)CLOCKS_PER_SEC;

	{
		unsigned nodeCount = Tree.nodeCount();

		if (nodeCount != Points.size())
			std::cout << "Inserted " << Points.size() << " points, but tree has only " << nodeCount << std::endl;
	}

	std::cout << "Tree with " << Points.size() << " points took " << runTime << " seconds to build." << std::endl;

	std::cout << "Tree height is " << Tree.TreeHeight() << ", there are " << Tree.nodeCount(true) << " nodes." << std::endl;

	performance(Tree, Points);

	Tree.rebuild();

	std::cout << "Rebuilt tree has " << Points.size() << " points, its height is " << Tree.TreeHeight() <<
				 ", and " << Tree.nodeCount(true) << " nodes." << std::endl;

	std::cout << std::endl;
	*/
	// Test tree built the efficient way

	std::cout << "Phase II\n";
	std::cout << "Building tree" << std::endl;

	startTime = clock();
	Tree.insert(Points);
	endTime = clock();

	runTime = (endTime - startTime) / static_cast<double>(CLOCKS_PER_SEC);
	
	if (Tree.nodeCount() != Points.size())
		std::cout << "Inserted " << Points.size() << " points, but tree has only " << Tree.nodeCount() << std::endl;

	std::cout << "Tree with " << Points.size() << " points took " << runTime << " seconds to build.\n";

	std::cout << "Tree height is " << Tree.TreeHeight() << ", there are " << Tree.nodeCount(true) << " nodes." << std::endl;

	performance(Tree);

	Tree.clear();
	/*
	std::cout << "Phase III" << std::endl;
	std::cout << "Building tree" << std::endl;

	for (size_t i = 0; i < Points.size(); i++)
		Tree.insert(Points[i]);

	unsigned nodes = Tree.nodeCount();

	if (nodes != Points.size())
		std::cout << "Err: inserted " << Points.size() << " points, but the tree has " << nodes << std::endl;
	else
		std::cout << "OK: Tree size matches number of inserted points." << std::endl;

	std::cout << "Erasing tree" << std::endl;

	std::random_shuffle(Points.begin(), Points.end());

	for (size_t i = 0; i < Points.size(); i++)
		Tree.erase(Points[i]);
	
	nodes = Tree.nodeCount(true);

	if (nodes == 0)
		std::cout << "OK: tree completely erased." << std::endl;
	else
		std::cout << "Err: erased tree has " << nodes << " nodes." << std::endl;
	*/
	// -----------------------------------------------------------------------------------------------------------------------

	return 0;
}
