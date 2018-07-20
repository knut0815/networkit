/*
 * Kadabra.h
 *
 * Created on: 18.07.2018
 * 		 Author: Eugenio Angriman, Alexander van der Grinten
 */

#ifndef KADABRA_H_
#define KADABRA_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

class Kadabra : public Algorithm {
public:
	Kadabra(const Graph &G, const count k = 0, const double delta = 0.1,
	        const double err = 0.01, const count seed = 42,
	        const count unionSample = 0, const count startFactor = 100);
	void run() override;
	std::vector<node> topkNodesList();
	std::vector<double> topkScoresList();

protected:
	const Graph &G;
	const count k;
	const double delta;
	const double err;
	const count n;
	const count seed;
	const count startFactor;

	count unionSample;
	count nPairs;

	std::vector<count> approx;
	std::vector<node> topkNodes;
	std::vector<double> topkScores;

	void init();
	std::vector<node> randomPath();
	void backtrackPath(const node u, const node v, const node start,
	                   std::vector<node> &path, std::vector<count> &nPaths,
	                   Graph &pred);
	void removeSomeEdges(Graph &pred, const std::vector<node> &vertices,
	                     const count length);

	count getDegree(const Graph &graph, node z, bool useDegreeIn);
};

inline std::vector<node> Kadabra::topkNodesList() { return topkNodes; }

inline std::vector<double> Kadabra::topkScoresList() { return topkScores; }

} // namespace NetworKit

#endif /* ifndef KADABRA_H_ */
