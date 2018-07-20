/*
 * Kadabra.cpp
 *
 * Created on: 18.07.2018
 * 		 Author: Eugenio Angriman, Alexander van der Grinten
 */

#include <cmath>
#include <omp.h>

#include "../auxiliary/PrioQueue.h"
#include "../auxiliary/Random.h"
#include "../distance/APSP.h"
#include "../distance/Diameter.h"
#include "Kadabra.h"

namespace NetworKit {
Kadabra::Kadabra(const Graph &G, const count k, const double delta,
                 const double err, const count seed, const count unionSample,
                 const count startFactor)
    : G(G), k(k), delta(delta), err(err), n(G.upperNodeIdBound()), seed(seed),
      startFactor(startFactor), unionSample(unionSample) {
	if (k > n) {
		throw std::runtime_error(
		    "k is higher than the number of nodes of the input graph!");
	}

	if (delta >= 1 || delta <= 0) {
		throw std::runtime_error(
		    "Delta should be greater than 0 and smaller than 1.");
	}

	if (err >= 1 || err <= 0) {
		throw std::runtime_error(
		    "The error should be greater than 0 and smaller than 1.");
	}
}

void Kadabra::init() {
	nPairs = 0;
	approx.assign(n, 0);
}

void Kadabra::run() {
	init();

	// TODO: check if exact diameter is required or if an upper bound is correct.
	// If so, which is the maximum relative error?
	Diameter diam(G, estimatedRange, 0);
	diam.run();
	// Getting diameter upper bound
	int32_t diameter = diam.getDiameter().second;
	INFO("Diameter upper bound is ", diameter);
	INFO("Diameter lower bound is ", diam.getDiameter().first);

	const double omega =
	    0.5 / err / err * (std::log2(diameter - 1) + 1 + std::log(0.5 / delta));

	const count tau = omega / startFactor;

	INFO("Omega is ", omega);
	INFO("Tau is ", tau);

	if (unionSample == 0) {
		unionSample =
		    std::min(n, (count)std::max((2 * std::sqrt(G.numberOfEdges()) /
		                                 omp_get_max_threads()),
		                                k + 20.));
	}

	Aux::PrioQueue<count, node> top(n);

#pragma omp parallel num_threads(1)
	{
		while (nPairs < tau) {
			auto path = randomPath();
			INFO(path);
#pragma omp critical
			{
				++nPairs;
				for (node u : path) {
					++approx[u];
					top.changeKey(approx[u], u);
				}
			}
			break; // continue here
		}
	}
}

std::vector<node> Kadabra::randomPath() {
	Graph pred(n, false, G.isDirected());

	node u = G.randomNode(), v = G.randomNode();
	while (u == v) {
		v = G.randomNode();
	}

	INFO("Random nodes: ", u, ", ", v);
	APSP apsp(G);
	apsp.run();
	INFO("Their distance is ", apsp.getDistance(u, v));

	std::vector<node> q(n, 0);
	count endQ = 2;
	q[0] = u;
	q[1] = v;

	std::vector<short> ballInd(n, 0);

	ballInd[u] = 1;
	ballInd[v] = 2;

	std::vector<count> dist(n, std::numeric_limits<count>::infinity());
	dist[u] = dist[v] = 0;

	std::vector<count> nPaths(n, 0);
	nPaths[u] = nPaths[v] = 1;

	std::vector<std::pair<node, node>> spEdges;

	node x, randomEdge;
	bool hasToStop = false;
	bool useDegreeIn;
	count startU = 0, startV = 1, endU = 1, endV = 2, startCur, endCur,
	      *newEndCur;
	count sumDegsU = G.degree(u), sumDegsV = G.degreeIn(v), *sumDegsCur;
	count totWeight = 0, curEdge = 0;

	auto procNeighbor = [&](node y) {
		if (ballInd[y] == 0) {
			(*sumDegsCur) += getDegree(G, y, useDegreeIn);
			nPaths[y] = nPaths[x];
			ballInd[y] = ballInd[x];
			q[endQ++] = y;
			(*newEndCur)++;
			pred.addEdge(y, x);
			dist[y] = dist[x] + 1;
		} else if (ballInd[x] != ballInd[y]) {
			hasToStop = true;
			spEdges.push_back(std::make_pair(x, y));
		} else if (dist[y] == dist[x] + 1) {
			nPaths[y] += nPaths[x];
			pred.addEdge(y, x);
		}
	};

	while (!hasToStop) {
		if (sumDegsU <= sumDegsV) {
			startCur = startU;
			endCur = endU;
			startU = endQ;
			newEndCur = &endU;
			endU = endQ;
			sumDegsU = 0;
			sumDegsCur = &sumDegsU;
			useDegreeIn = false;
		} else {
			startCur = startV;
			endCur = endV;
			startV = endQ;
			newEndCur = &endV;
			endV = endQ;
			sumDegsV = 0;
			sumDegsCur = &sumDegsV;
			useDegreeIn = true;
		}

		while (startCur < endCur) {
			x = q[startCur++];

			if (useDegreeIn) {
				G.forInNeighborsOf(x, [&](node y) { procNeighbor(y); });
			} else {
				G.forNeighborsOf(x, [&](node y) { procNeighbor(y); });
			}
		}

		if (*sumDegsCur == 0) {
			hasToStop = true;
		}
	}

	if (spEdges.size() == 0) {
		removeSomeEdges(pred, q, endQ);
		return std::vector<node>();
	}

	for (auto p : spEdges) {
		totWeight += nPaths[p.first] * nPaths[p.second];
	}

	randomEdge = Aux::Random::integer(totWeight);
	std::vector<node> path;

	for (auto p : spEdges) {
		curEdge += nPaths[p.first] * nPaths[p.second];
		if (curEdge > randomEdge) {
			backtrackPath(u, v, p.first, path, nPaths, pred);
			backtrackPath(u, v, p.second, path, nPaths, pred);
			break;
		}
	}

	return path;
}

void Kadabra::backtrackPath(const node u, const node v, const node start,
                            std::vector<node> &path, std::vector<count> &nPaths,
                            Graph &pred) {

	if (start == u || start == v) {
		return;
	}

	count totWeight = nPaths[start];
	node randomPred, curPred = 0;
	node w = 0;

	path.push_back(start);
	randomPred = Aux::Random::integer(totWeight);

	for (node t : pred.neighbors(start)) {
		w = t;
		curPred += nPaths[v];
		if (curPred > randomPred) {
			break;
		}
	}

	if (w != u && w != v) {
		backtrackPath(u, v, w, path, nPaths, pred);
	}
}

void Kadabra::removeSomeEdges(Graph &pred, const std::vector<node> &vertices,
                              const count length) {
	node u;
	for (count i = 0; i < length; ++i) {
		u = vertices[i];
		pred.forEdgesOf(u, [&](node v) { pred.removeEdge(u, v); });
		pred.forInEdgesOf(u, [&](node v) { pred.removeEdge(u, v); });
	}
}

count Kadabra::getDegree(const Graph &graph, node z, bool useDegreeIn) {
	return useDegreeIn ? graph.degreeIn(z) : graph.degree(z);
}

} // namespace NetworKit
