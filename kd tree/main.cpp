// kd tree.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "kd tree.h"

int _tmain(int argc, _TCHAR* argv[])
{
	kd_tree<3> Tree;
	kd_tree<3>::kd_point Point = { 0.5f, 0.5f, 0.5f },
						 MinCorner = { 0.0f, 0.0f, 0.0f },
						 MaxCorner = { 2.0f, 2.0f, 2.0f };
	std::vector<kd_tree<3>::kd_point> Points, nearPoints;

	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 0.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 0.0f, 1.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 1.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 1.0f, 1.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 1.0f, 0.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 1.0f, 0.0f, 1.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 1.0f, 1.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 1.0f, 1.0f, 1.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 0.0f, 2.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 2.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 0.0f, 2.0f, 2.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 2.0f, 0.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 2.0f, 0.0f, 2.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 2.0f, 2.0f, 0.0f }));
	Points.push_back(kd_tree<3>::kd_point({ 2.0f, 2.0f, 2.0f }));

	Tree.insert(Points);

	Tree.PrintTree();

	std::cout << std::endl;

	std::cout << "Tree height is " << Tree.TreeHeight() << '\n';
	std::cout << "There are " << Tree.nodeCount() << " points in the tree.\n";

	std::cout << std::endl;
	std::cout << "Points in cloud are:\n";

	for (kd_tree<3>::kd_point Point : Tree)
	{
		std::cout << Point[0] << " " << Point[1] << " " << Point[2] << '\n';
	}
	std::cout << std::endl;

	return 0;
}
