#include <array>
#include <vector>
#include <limits>
#include <memory>
#include <optional>
#include <iterator>
#include <iostream>
#include <algorithm>

/***************************************************/
/*												   */
/* 			   (C) Uri Raz, 2015-2016			   */
/*			   uri.raz@private.org.il			   */
/*												   */
/* Shared for the greater glory of the C.P.O.N.O.D */
/*												   */
/***************************************************/

// Revision history
//  2015-11-12 First stable version, reviewed & optimized.
//  2015-11-24 Replaced KNearestNeighbors's use of unordered_set with linear search
//  2015-12-24 Cleanups following user feedback
//  2016-01-15 Add point insertion, point deletion, and tree rebuilding.
//  2016-10-12 Refactoring & improvements
//  2017-11-04 Cleanups
//  2017-12-14 Cleanups

// ---------------------------------------------------------------------------------------------

template <unsigned int K>
unsigned diffIndex(const std::array<float, K> &P, const std::array<float, K> &Q, unsigned s = 0)
{
	for (unsigned i = 0; i < K; ++i, ++s %= K)
		if (P[s] != Q[s])
			return s;

	return K;
}

// -------------------------------------------------------

template <unsigned int K>
class kd_box
{
public:
	using kd_point = std::array<float, K>;
	kd_point first = { 0 }, second = { 0 };
	
	kd_box(const kd_point &Point1, const kd_point &Point2)
	{
		for (unsigned i = 0; i < K; i++)
		{
			first[i]  = std::min(Point1[i], Point2[i]);
			second[i] = std::max(Point1[i], Point2[i]);
		}
	}

	kd_box(const kd_box &Box1, const kd_box &Box2)
	{
		for (unsigned i = 0; i < K; i++)
		{
			first[i]  = std::min(Box1.first[i], Box2.first[i]);
			second[i] = std::max(Box1.second[i], Box2.second[i]);
		}
	}
	
	kd_box(const kd_point &Point, const float distance)
	{
		for (unsigned i = 0; i < K; i++)
		{
			first[i]  = Point[i] - distance;
			second[i] = Point[i] + distance;
		}
	}

	kd_box(const typename std::vector<kd_point>::iterator &iterBegin,
		   const typename std::vector<kd_point>::iterator &iterEnd)
	{
		if (iterBegin != iterEnd)
		{
			first = second = *iterBegin;

			for (auto currPoint = std::next(iterBegin); currPoint != iterEnd; currPoint++)
			{
				for (unsigned i = 0; i < K; i++)
				{
					first[i]  = std::min(first[i], (*currPoint)[i]);
					second[i] = std::max(second[i], (*currPoint)[i]);
				}
			}
		}
		else
		{
			first.fill(std::numeric_limits<float>::quiet_NaN());
			second.fill(std::numeric_limits<float>::quiet_NaN());
		}
	}

	void updateBox(const kd_point &Point)
	{
		for (unsigned i = 0; i < K; i++)
		{
			first[i]  = std::min(first[i], Point[i]);
			second[i] = std::max(second[i], Point[i]);
		}
	}
};

template <unsigned int K>
inline bool operator==(const kd_box<K>& lhs, const kd_box<K>& rhs)
{
	return (lhs.first  == rhs.first &&
			lhs.second == rhs.second);
}

template <unsigned int K>
unsigned diffIndex(const kd_box<K> &box, unsigned s = 0)
{
	return diffIndex(box.first, box.second, s);
}

// ---------------------------------------------------------------------------------------------

template <unsigned int K>
class kd_tree
{
public:

	using kd_point = std::array<float, K>;

private:

	static bool isValid(const kd_point &kdPoint)
	{
		for (auto i : kdPoint)
			if (!std::isfinite(i))
				return false;

		return true;
	}

protected:

	class kd_node
	{
	public:
		virtual bool isInternal() const = 0;
		virtual kd_box<K> boundingBox() const = 0;
		virtual unsigned TreeHeight() const = 0;

		virtual void pointsInBox(const kd_box<K> &searchBox, std::vector<kd_point> &Points) const = 0;

		virtual void nearestNeighbor(const kd_point &srcPoint, kd_point &nearPoint,
									 float &minDistance, kd_box<K> &minRegion) const = 0;

		virtual void KNearestNeighbors(const kd_point &srcPoint, const unsigned k,
									   std::vector<kd_point> &nearPoints, std::vector<float> &nearDistances,
									   kd_box<K> &minRegion) const = 0;

		virtual unsigned nodeCount(bool withInternalNodes) const = 0;
	};

	// ------------------------------------------

	class kd_internal_node final : public kd_node
	{
	private:
		float m_splitVal;
		unsigned m_axis;
		kd_box<K> m_boundingBox;

	public:
		std::shared_ptr<kd_node> m_Left, m_Right;

		kd_internal_node(const float splitVal, const unsigned axis, const kd_box<K> &boundingBox,
						 const std::shared_ptr<kd_node> &Left, const std::shared_ptr<kd_node> &Right) :
									m_splitVal(splitVal), m_axis(axis), m_boundingBox(boundingBox), m_Left(Left), m_Right(Right) {}

		void updateBox(const kd_point &Point) { m_boundingBox.updateBox(Point);  }

		unsigned splitAxis() const { return m_axis; }

		float splitVal() const { return m_splitVal; }

		virtual bool isInternal() const override  { return true; }

		virtual kd_box<K> boundingBox() const override  { return m_boundingBox; }

		virtual void pointsInBox(const kd_box<K> &searchBox, std::vector<kd_point> &Points) const override 
		{
			if (regionCrossesRegion(searchBox, m_Left->boundingBox()))
				m_Left->pointsInBox(searchBox, Points);

			if (regionCrossesRegion(searchBox, m_Right->boundingBox()))
				m_Right->pointsInBox(searchBox, Points);
		}

		virtual void KNearestNeighbors(const kd_point &srcPoint, const unsigned k,
									   std::vector<kd_point> &nearPoints, std::vector<float> &nearDistances,
									   kd_box<K> &minRegion) const override 
		{
			if (regionCrossesRegion(m_Left->boundingBox(), minRegion))
				m_Left->KNearestNeighbors(srcPoint, k, nearPoints, nearDistances, minRegion);

			if (regionCrossesRegion(m_Right->boundingBox(), minRegion))
				m_Right->KNearestNeighbors(srcPoint, k, nearPoints, nearDistances, minRegion);
		}

		virtual void nearestNeighbor(const kd_point &srcPoint, kd_point &nearPoint,
									 float &minDistance, kd_box<K> &minRegion) const override 
		{
			if (regionCrossesRegion(m_Left->boundingBox(), minRegion))
				m_Left->nearestNeighbor(srcPoint, nearPoint, minDistance, minRegion);

			if (regionCrossesRegion(m_Right->boundingBox(), minRegion))
				m_Right->nearestNeighbor(srcPoint, nearPoint, minDistance, minRegion);
		}

		virtual unsigned TreeHeight() const override 
		{
			return 1 + std::max(m_Left->TreeHeight(), m_Right->TreeHeight());
		}

		virtual unsigned nodeCount(bool withInternalNodes) const override 
		{
			return (withInternalNodes ? 1 : 0) + m_Left->nodeCount(withInternalNodes) + m_Right->nodeCount(withInternalNodes);
		}
	};		// kd_internal_node

	// --------------------------------------

	class kd_leaf_node final : public kd_node
	{
	private:
		kd_point m_pointCoords;

	public:
		std::weak_ptr<kd_leaf_node> m_Prev, m_Next;

		kd_leaf_node(const kd_point &Point) : m_pointCoords(Point) { }

		kd_leaf_node(const kd_point &Point, std::weak_ptr<kd_leaf_node> Prev, std::weak_ptr<kd_leaf_node> Next) :
			m_pointCoords(Point), m_Prev(Prev), m_Next(Next) {}

		const kd_point pointCoords() const { return m_pointCoords; }

		virtual bool isInternal() const override  { return false; }

		virtual unsigned TreeHeight() const override  { return 1; }

		virtual unsigned nodeCount(bool) const override  { return 1; }

		virtual kd_box<K> boundingBox() const override  { return { m_pointCoords, m_pointCoords }; }

		virtual void pointsInBox(const kd_box<K>&, std::vector<kd_point> &Points) const override
		{
			Points.push_back(m_pointCoords);
		}

		// ----------------------------------------------------------------------------------

		virtual void nearestNeighbor(const kd_point &srcPoint, kd_point &nearPoint,
									 float &minDistance, kd_box<K> &minRegion) const override
		{
			float currDistance = Distance(srcPoint, m_pointCoords);

			if (currDistance < minDistance)
			{
				nearPoint   = m_pointCoords;
				minDistance = currDistance;
				minRegion   = kd_box<K>(srcPoint, minDistance);
			}
		}

		// -------------------------------------------------------------------------------------------------

		virtual void KNearestNeighbors(const kd_point &srcPoint, const unsigned k,
									   std::vector<kd_point> &nearPoints, std::vector<float> &nearDistances,
									   kd_box<K> &minRegion) const override
		{
			float    currDistance = Distance(srcPoint, m_pointCoords);
			unsigned i = k - 1;

			if (currDistance > nearDistances[i]) return;

			for (const kd_point& elem : nearPoints)
				if (elem == m_pointCoords) return;

			nearPoints[i] = m_pointCoords;
			nearDistances[i] = currDistance;

			while (i > 0 && nearDistances[i - 1] > nearDistances[i])
			{
				std::swap(nearPoints[i - 1], nearPoints[i]);
				std::swap(nearDistances[i - 1], nearDistances[i]);
				i--;
			}

			minRegion = kd_box<K>(srcPoint, nearDistances[k - 1]);
		}
	};		// kd_leaf_node

	// ---------------------------------------------------------------------

	// The routine has a desired side effect of sorting Points
	float NthCoordMedian(typename std::vector<kd_point>::iterator iterBegin,
						 typename std::vector<kd_point>::iterator iterEnd,
						 const unsigned num)
	{
		std::sort(iterBegin, iterEnd, [num](const kd_point &A, const kd_point &B) { return A[num] < B[num]; });

		const auto numPoints = iterEnd - iterBegin;

		float Median = (*(iterBegin + numPoints / 2))[num];

		if (numPoints % 2 == 0)
			Median = (Median + (*(iterBegin + numPoints / 2 - 1))[num]) / 2.0f;

		if (Median == (*(iterEnd - 1))[num] && Median != (*iterBegin)[num])
		{
			auto medianIter = iterBegin + numPoints / 2;

			while (Median == (*(--medianIter))[num]);

			Median = (Median + (*medianIter)[num]) / 2.0f;
		}

		return Median;
	}

	// ------------------------------------------------------------------------------------

	std::shared_ptr<kd_node> CreateTree(typename std::vector<kd_point>::iterator iterBegin,
										typename std::vector<kd_point>::iterator iterEnd,
										std::shared_ptr<kd_leaf_node> &lastLeaf,
										unsigned Depth = 0)
	{
		if (iterEnd - iterBegin == 1)
		{
			std::shared_ptr<kd_leaf_node> retNode(std::make_shared<kd_leaf_node>(*iterBegin));

			if (lastLeaf)
			{
				lastLeaf->m_Next = retNode;
				retNode->m_Prev = lastLeaf;
			}
			else
				m_firstLeaf = retNode;

			lastLeaf = retNode;

			return retNode;
		}
		else if (iterEnd - iterBegin == 2)
		{
			kd_point point0 = *iterBegin, point1 = *(std::next(iterBegin));
			unsigned splitAxis = diffIndex(point0, point1, Depth%K);

			if (point0[splitAxis] > point1[splitAxis])
				std::swap(point0, point1);

			std::shared_ptr<kd_leaf_node> Left(std::make_shared<kd_leaf_node>(point0)), Right(std::make_shared<kd_leaf_node>(point1));

			if (lastLeaf)
			{
				lastLeaf->m_Next = Left;
				Left->m_Prev = lastLeaf;
			}
			else
				m_firstLeaf = Left;

			Left->m_Next = Right;
			Right->m_Prev = Left;

			lastLeaf = Right;

			std::shared_ptr<kd_internal_node> retNode(std::make_shared<kd_internal_node>(
												(point0[splitAxis] + point1[splitAxis]) / 2.0f, splitAxis,kd_box<K>(point0, point1), Left, Right));

			return retNode;
		}
		else
		{
			kd_box<K> boundingBox(iterBegin, iterEnd);

			const unsigned  splitAxis = diffIndex(boundingBox, Depth%K);
			const float     Median    = NthCoordMedian(iterBegin, iterEnd, splitAxis);
			const size_t    numPoints = iterEnd - iterBegin;

			auto lastMedianLoc = iterBegin + numPoints / 2;

			if ((*lastMedianLoc)[splitAxis] != (*std::prev(iterEnd))[splitAxis])
				while ((*(++lastMedianLoc))[splitAxis] == Median);
			else
				while ((*std::prev(lastMedianLoc))[splitAxis] == Median)
					lastMedianLoc--;

			std::shared_ptr<kd_node> Left, Right;
			
			Left  = CreateTree(iterBegin, lastMedianLoc, lastLeaf, Depth + 1);
			Right = CreateTree(lastMedianLoc, iterEnd,   lastLeaf, Depth + 1);

			std::shared_ptr<kd_internal_node> retNode(std::make_shared<kd_internal_node>(Median, splitAxis, boundingBox, Left, Right));

			return retNode;
		}
	}

	// ------------------------------------------------------------------------------------

	std::shared_ptr<kd_leaf_node> ApproxNearestNeighborNode(const kd_point &srcPoint) const
	{
		std::shared_ptr<kd_node> Node(m_Root);
		std::shared_ptr<kd_internal_node> iNode;

		while (Node->isInternal())
		{
			iNode = std::static_pointer_cast<kd_internal_node>(Node);

			Node = (srcPoint[iNode->splitAxis()] <= iNode->splitVal()) ? iNode->m_Left : iNode->m_Right;
		}

		return std::static_pointer_cast<kd_leaf_node>(Node);
	}

	// ----------------------------------------------------------------

	kd_point ApproxNearestNeighborPoint(const kd_point &srcPoint) const
	{
		return ApproxNearestNeighborNode(srcPoint)->pointCoords();
	}

	// ----------------------------------------------------------------

	float ApproxNearestNeighborDistance(const kd_point &srcPoint) const
	{
		return Distance(srcPoint, ApproxNearestNeighborPoint(srcPoint));
	}

	// -------------------------------------

	std::shared_ptr<kd_node> m_Root;
	std::weak_ptr<kd_leaf_node> m_firstLeaf;

public:

	kd_tree() { }

	kd_tree(std::vector<kd_point> &Points) { insert(Points); }

	void clear() { m_Root.reset(); }

	kd_tree(const kd_tree &obj) = delete;
	bool operator=(const kd_tree<K> &rhs) = delete;
	bool operator==(const kd_tree<K> rhs) = delete;

	std::optional<kd_box<K>> boundingBox() const
	{
		if (m_Root)
			return m_Root->boundingBox();
		else
			return std::nullopt;
	}

	friend void swap(kd_tree& a, kd_tree& b)
	{
		using std::swap;

		swap(a.m_Root, b.m_Root);
		swap(a.m_firstLeaf, b.m_firstLeaf);
	}

	static bool pointIsInRegion(const kd_point &Point, const kd_box<K> &Region);
	static bool regionCrossesRegion(const kd_box<K> &Region1, const kd_box<K> &Region2);

	static float Distance(const kd_point &P, const kd_point &Q);
	static float DistanceSq(const kd_point &P, const kd_point &Q);

	// ---------------------------------------------------------------------------------------------------------
	//
	// Iterator implementation
	//
	class const_iterator : public std::iterator<std::forward_iterator_tag, kd_point, void, kd_point*, kd_point&>
	{
	public:

		const_iterator() = default;
		const_iterator(const std::shared_ptr<kd_leaf_node> &node) : nodePtr(node) {}

		bool operator==(const const_iterator &rhs) { return this->nodePtr.lock() == rhs.nodePtr.lock(); }
		bool operator!=(const const_iterator &rhs) { return this->nodePtr.lock() != rhs.nodePtr.lock(); }

		kd_point operator*() { return this->nodePtr.lock()->pointCoords(); }

		const_iterator& operator++()
		{
			this->nodePtr = this->nodePtr.lock()->m_Next.lock();

			return *this;
		}

		const_iterator operator++(int)
		{
			const_iterator tmp(*this);
			++*this;
			return tmp;
		}

	private:
		std::weak_ptr<kd_leaf_node> nodePtr;
	};

	const_iterator begin()
	{
		return m_Root ? m_firstLeaf.lock() : end();
	}

	const_iterator end()
	{
		const_iterator retVal;

		return retVal;
	}
	//
	// End of Iterator implementation
	//
	// ------------------------------

	// ---------------------------------------

	void insert(std::vector<kd_point> &Points)
	{
		this->clear();

		for (signed i = Points.size() - 1; i >= 0; i--)
			if (!isValid(Points[i]))
				Points.erase(Points.begin() + i);

		if (Points.size() > 0)
		{
			sort(Points.begin(), Points.end());
			auto it = unique(Points.begin(), Points.end());
			Points.resize(distance(Points.begin(), it));
			Points.shrink_to_fit();

			std::shared_ptr<kd_leaf_node> dummy;

			m_Root = CreateTree(Points.begin(), Points.end(), dummy);
		}
	}

	// ---------------------------------------------------------------------------------------------

	bool insert(const kd_point &Point)
	{
		if (!isValid(Point)) return false;

		if (!m_Root)
		{
			m_Root = std::make_shared<kd_leaf_node>(Point);

			m_firstLeaf = std::static_pointer_cast<kd_leaf_node>(m_Root);
		}
		else if (!m_Root->isInternal())
		{
			std::shared_ptr<kd_leaf_node> rootNode = std::static_pointer_cast<kd_leaf_node>(m_Root);

			kd_point rootPoint = rootNode->pointCoords();

			if (Point == rootPoint) return false;

			unsigned splitAxis = diffIndex(rootPoint, Point);
			kd_point prmPoint = Point;

			if (Point[splitAxis] > rootPoint[splitAxis])
				std::swap(prmPoint, rootPoint);

			std::shared_ptr<kd_leaf_node> Left(std::make_shared<kd_leaf_node>(prmPoint)), Right(std::make_shared<kd_leaf_node>(rootPoint));

			Left->m_Next = Right;
			Right->m_Prev = Left;

			m_Root = std::make_shared<kd_internal_node>
				((rootPoint[splitAxis] + prmPoint[splitAxis]) / 2.0f, splitAxis, kd_box<K>(rootPoint, prmPoint), Left, Right);

			m_firstLeaf = Left;
		}
		else
		{
			std::shared_ptr<kd_node> currNode(m_Root);
			std::shared_ptr<kd_internal_node> prevNode;

			while (currNode->isInternal())
			{
				prevNode = std::static_pointer_cast<kd_internal_node>(currNode);

				prevNode->updateBox(Point);

				currNode = Point[prevNode->splitAxis()] <= prevNode->splitVal() ? prevNode->m_Left : prevNode->m_Right;
			}

			const bool prvSide = (Point[prevNode->splitAxis()] <= prevNode->splitVal());

			std::shared_ptr<kd_leaf_node> brother =
								std::static_pointer_cast<kd_leaf_node>(prvSide ? prevNode->m_Left : prevNode->m_Right);

			kd_point brotherPoint = brother->pointCoords();

			if (Point == brotherPoint) return false;

			unsigned splitAxis = diffIndex(Point, brotherPoint);

			std::shared_ptr<kd_leaf_node> pointNode(std::make_shared<kd_leaf_node>(Point));
			std::shared_ptr<kd_internal_node> iNode;

			if (Point[splitAxis] <= brotherPoint[splitAxis])
			{
				iNode = std::make_shared<kd_internal_node>((Point[splitAxis] + brotherPoint[splitAxis]) / 2.0f, splitAxis,
														   kd_box<K>(Point, brotherPoint), pointNode, brother);

				pointNode->m_Next = brother;
				pointNode->m_Prev = brother->m_Prev;

				brother->m_Prev = pointNode;

				std::shared_ptr<kd_leaf_node> pPrev = pointNode->m_Prev.lock();

				if (pPrev)
					pPrev->m_Next = pointNode;

				if (m_firstLeaf.lock() == brother)
					m_firstLeaf = pointNode;
			}
			else
			{
				iNode = std::make_shared<kd_internal_node>((Point[splitAxis] + brotherPoint[splitAxis]) / 2.0f, splitAxis,
					kd_box<K>(Point, brotherPoint), brother, pointNode);

				pointNode->m_Next = brother->m_Next;
				pointNode->m_Prev = brother;

				brother->m_Next = pointNode;

				std::shared_ptr<kd_leaf_node> pNext = pointNode->m_Next.lock();

				if (pNext)
					pNext->m_Prev = pointNode;
			}

			(prvSide ? prevNode->m_Left : prevNode->m_Right) = iNode;
		}

		return true;
	}

	// ---------------------------------------------------------------------------------------------

	bool erase(const kd_point &Point)
	{
		if (m_Root == nullptr || !isValid(Point)) return false;

		if (!m_Root->isInternal())
		{
			std::shared_ptr<kd_leaf_node> rootNode = std::static_pointer_cast<kd_leaf_node>(m_Root);

			if (Point != rootNode->pointCoords()) return false;

			clear();
		}
		else
		{
			std::shared_ptr<kd_node> Node(m_Root);
			std::shared_ptr<kd_internal_node> p1Node, p2Node;

			while (Node->isInternal())
			{
				p2Node = p1Node;
				p1Node = std::static_pointer_cast<kd_internal_node>(Node);

				Node = Point[p1Node->splitAxis()] <= p1Node->splitVal() ? p1Node->m_Left : p1Node->m_Right;
			}

			std::shared_ptr<kd_leaf_node> lNode = std::static_pointer_cast<kd_leaf_node>(Node);

			if (Point != lNode->pointCoords()) return false;

			// Update the iterator's linked list

			std::shared_ptr<kd_leaf_node> pNode = lNode->m_Prev.lock(), nNode = lNode->m_Next.lock();

			if (pNode == nullptr)
				m_firstLeaf = nNode;
			else if (nNode)
			{
				pNode->m_Next = lNode->m_Next;
				nNode->m_Prev = lNode->m_Prev;
			}

			// Remove point from tree

			if (p2Node)
				( (Point[p2Node->splitAxis()] <= p2Node->splitVal() ) ? p2Node->m_Left  : p2Node->m_Right ) =
				   Point[p1Node->splitAxis()] <= p1Node->splitVal()   ? p1Node->m_Right : p1Node->m_Left;
			else
				m_Root = Point[p1Node->splitAxis()] <= p1Node->splitVal() ?
							std::static_pointer_cast<kd_internal_node>(m_Root)->m_Right :
							std::static_pointer_cast<kd_internal_node>(m_Root)->m_Left;
		}

		return true;
	}

	// ----------------------------------

	void rebuild()
	{
		unsigned leafNodes = nodeCount();

		if (leafNodes < 3) return;

		std::vector<kd_point> points;

		points.reserve(leafNodes);
		points.shrink_to_fit();

		for (const kd_point& Point : *this)
			points.push_back(Point);

		insert(points);
	}

	// --------------------------------------------------------------------------------

	void pointsInBox(const kd_box<K> &searchBox, std::vector<kd_point> &Points) const
	{
		Points.clear();

		if (m_Root && regionCrossesRegion(searchBox, m_Root->boundingBox()))
			m_Root->pointsInBox(searchBox, Points);

		Points.shrink_to_fit();
	}

	// ----------------------------------------------------------------------

	std::optional<kd_point> nearestNeighbor(const kd_point &srcPoint) const
	{
		if (!m_Root) return std::nullopt;

		kd_point nearPoint = ApproxNearestNeighborPoint(srcPoint);
		float minDistance = Distance(srcPoint, nearPoint);
		
		auto minBox = kd_box<K>(srcPoint, minDistance);
		m_Root->nearestNeighbor(srcPoint, nearPoint, minDistance, minBox);

		return nearPoint;
	}

	// --------------------------------------------------------------------------------------------------------

	std::vector<kd_point> KNearestNeighbors(const kd_point &srcPoint, const unsigned k) const
	{
		std::vector<kd_point> nearPoints;

		if (!m_Root) nearPoints;

		std::shared_ptr<kd_leaf_node> nNode = ApproxNearestNeighborNode(srcPoint),
									  pNode = nNode->m_Prev.lock();

		nearPoints.reserve(k);
		nearPoints.push_back(nNode->pointCoords());
		nNode = nNode->m_Next.lock();

		while (nearPoints.size() < k && (nNode || pNode))
		{
			if (nNode)
				nearPoints.push_back(nNode->pointCoords()),
				nNode = nNode->m_Next.lock();

			if (pNode && nearPoints.size() < k)
				nearPoints.push_back(pNode->pointCoords()),
				pNode = pNode->m_Prev.lock();
		}

		sort(nearPoints.begin(), nearPoints.end(),
			 [srcPoint](const kd_point &A, const kd_point &B) {return DistanceSq(srcPoint, A) < DistanceSq(srcPoint, B); });

		if (nearPoints.size() == k)
		{
			std::vector<float> minDistances;
			minDistances.reserve(k);
			minDistances.resize(k);

			for (unsigned i = 0; i < k; i++)
				minDistances[i] = Distance(srcPoint, nearPoints[i]);

			auto minBox = kd_box<K>(srcPoint, minDistances[k - 1]);
			m_Root->KNearestNeighbors(srcPoint, k, nearPoints, minDistances, minBox);
		}

		return nearPoints;
	}

	// ---------------------------------------------------------------------------

	void PrintTree(const std::shared_ptr<kd_node> &node, unsigned depth = 0) const
	{
		for (unsigned i = 0; i < depth; i++) std::cout << ' ';

		if (node == nullptr)
			std::cout << "null" << std::endl;
		else
		{
			if (node->isInternal())
			{
				std::shared_ptr<kd_internal_node> iNode = std::static_pointer_cast<kd_internal_node>(node);

				std::cout << "Split val is " << iNode->splitVal() << " for axis #" << iNode->splitAxis() << '\n';

				PrintTree(iNode->m_Left,  depth + 1);
				PrintTree(iNode->m_Right, depth + 1);
			}
			else
			{
				std::shared_ptr<kd_leaf_node> lNode = std::static_pointer_cast<kd_leaf_node>(node);

				kd_point point = lNode->pointCoords();

				std::cout << "Point is (";

				for (const float &val : point)
					std::cout << val << " ";

				std::cout << ")" << '\n';
			}
		}
	}

	// -----------------------------------------------------

	unsigned nodeCount(bool withInternalNodes = false) const
	{
		return m_Root ? m_Root->nodeCount(withInternalNodes) : 0;
	}

	// -----------------------------------------------------------

	unsigned TreeHeight() const
	{
		return this->m_Root ? m_Root->TreeHeight() : 0;
	}

	// ------------------------------------------

	void PrintTree() const { PrintTree(m_Root); }
};

// -------------------------------------------------------------

template <unsigned int K>
float kd_tree<K>::Distance(const kd_point &P, const kd_point &Q)
{
	return std::sqrtf(DistanceSq(P, Q));
}

// ---------------------------------------------------------------

template <unsigned K>
float kd_tree<K>::DistanceSq(const kd_point &P, const kd_point &Q)
{
	float Sum = 0;

	for (unsigned i = 0; i < K; i++)
		Sum += (P[i] - Q[i]) * (P[i] - Q[i]);

	return Sum;
}

// -----------------------------------------------------------------------------

template <unsigned K>
bool kd_tree<K>::pointIsInRegion(const kd_point &Point, const kd_box<K> &Region)
{
	for (unsigned i = 0; i < K; i++)
		if (!(Region.first[i] <= Point[i] && Point[i] <= Region.second[i]))
			return false;

	return true;
}

// --------------------------------------------------------------------------------------

template <unsigned K>
bool kd_tree<K>::regionCrossesRegion(const kd_box<K> &Region1, const kd_box<K> &Region2)
{
	for (unsigned i = 0; i < K; i++)
		if (Region1.first[i] > Region2.second[i] || Region1.second[i] < Region2.first[i])
			return false;

	return true;
}
