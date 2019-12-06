#include <iostream>
#include <limits.h> 
#include <string.h>
#include <queue>
#include <array>
#include <random>
#include "boost/random/piecewise_linear_distribution.hpp"
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <fstream>
#include <vector>
#include <string>
#include <stack>
#include "..\..\..\gnuplot\gnuplot-iostream\gnuplot-i.h"
#include <complex>
#include <cmath>
#include <conio.h>   //for getch(), needed in wait_for_key()
#include <windows.h> //for Sleep()
#include "boost\tuple\tuple.hpp"
#include "boost\array.hpp"
#include "boost\range\adaptor\transformed.hpp"
#include "boost\range\irange.hpp"
#include "boost\bind.hpp"

using namespace std;
#ifdef USE_ARMA
#include <armadillo>
#endif

#ifdef USE_BLITZ
#include <blitz/array.h>
#endif


#define DBG 1   // set DBG 1 for debugging code and 0 for normal run
#define GNUPLOT_ENABLE_PTY

using namespace std;
// Number of vertices in given graph 
#define V 6 
#define TABSIZE 4

struct nodeStruct {
	string name;
	float duration;
	float cost;
	int es, ef, ls, lf, st;  // es : earliest start time , ef : earliest finish time
							 // ls : latest start time ,  lf : latest finish time
							 // st : slack time 
} node;


void wait_for_key()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
	cout << endl << "Press any key to continue..." << endl;

	FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
	_getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
	cout << endl << "Press ENTER to continue..." << endl;

	std::cin.clear();
	std::cin.ignore(std::cin.rdbuf()->in_avail());
	std::cin.get();
#endif
	return;
}
void test1(Gnuplot gp) {

	gp.reset_plot();
	gp.set_xlabel("Activity Duration");
	gp.set_ylabel("Activity Cost");
	gp.plotfile_xy("Data.txt", 1, 2, "OutputGraph");
	gp.showonscreen(); // window output
	wait_for_key();
}


void test(Gnuplot gp) {

	//
	// Equations
	//
	gp.reset_plot();
	cout << endl << endl << "*** various equations" << endl;
	gp.set_title("Slopes\\nNew Line");
	gp.set_xrange(1, 100);
	gp.set_yrange(1, 50);
	cout << "y = x" << endl;
	gp.plot_slope(1.0, 1.0, "y=x");



	gp.showonscreen(); // window output
	wait_for_key();

}

void plotprocesses() {

	Gnuplot gp;
	Gnuplot::set_GNUPlotPath("C:/gnuplot/bin");
	//test(gp);
	test1(gp);
	gp.showonscreen(); // window output
	wait_for_key();

}



// returns vector of n numbers for input
std::vector<int> ReadNumbers()
{
	std::vector<int> numbers;
	do
	{
		int input;
		if (std::cin >> input)
			numbers.push_back(input);
	} while (std::cin && std::cin.peek() != '\n');

	return numbers;
}

// utility for topological sorting of activity graph
void topologicalSortUtil(int v, vector<bool>& visited, stack<int>& Stack, vector< vector<int> >& adj)
{
	visited[v] = true;

	vector<int>::iterator i;
	for (i = adj[v].begin(); i != adj[v].end(); ++i)
		if (!visited[*i])
			topologicalSortUtil(*i, visited, Stack, adj);

	Stack.push(v);
}


/* Returns true if there is a path from source 's' to sink 't' in
  residual graph. Also fills parent[] to store the path */
int bfs(int rGraph[V][V], int s, int t, int parent[])
{
	// Create a visited array and mark all vertices as not visited 
	bool visited[V];
	memset(visited, 0, sizeof(visited));

	// Create a queue, enqueue source vertex and mark source vertex 
	// as visited 
	queue <int> q;
	q.push(s);
	visited[s] = true;
	parent[s] = -1;

	// Standard BFS Loop 
	while (!q.empty())
	{
		int u = q.front();
		q.pop();

		for (int v = 0; v < V; v++)
		{
			if (visited[v] == false && rGraph[u][v] > 0)
			{
				q.push(v);
				parent[v] = u;
				visited[v] = true;
			}
		}
	}

	// If we reached sink in BFS starting from source, then return 
	// true, else false 
	return (visited[t] == true);
}

// A DFS based function to find all reachable vertices from s.  The function 
// marks visited[i] as true if i is reachable from s.  The initial values in 
// visited[] must be false. We can also use BFS to find reachable vertices 
void dfs(int rGraph[V][V], int s, bool visited[])
{
	visited[s] = true;
	for (int i = 0; i < V; i++)
		if (rGraph[s][i] && !visited[i])
			dfs(rGraph, i, visited);
}

int fordFulkerson(int graph[V][V], int numOfNodes, int startNode, int endNode)
{
	int u, v;

	// Create a residual graph and fill the residual graph with 
	// given capacities in the original graph as residual capacities 
	// in residual graph 
	int rGraph[V][V]; // rGraph[i][j] indicates residual capacity of edge i-j 
	for (u = 0; u < V; u++)
		for (v = 0; v < V; v++)
			rGraph[u][v] = graph[u][v];

	int parent[V];  // This array is filled by BFS and to store path 
	int max_flow = 0;
	// Augment the flow while there is a path from source to sink 
	while (bfs(rGraph, startNode, endNode, parent))
	{
		// Find minimum residual capacity of the edhes along the 
		// path filled by BFS. Or we can say find the maximum flow 
		// through the path found. 
		int path_flow = INT_MAX;
		for (v = endNode; v != startNode; v = parent[v])
		{
			u = parent[v];
			path_flow = min(path_flow, rGraph[u][v]);
		}

		// update residual capacities of the edges and reverse edges 
		// along the path 
		for (v = endNode; v != startNode; v = parent[v])
		{
			u = parent[v];
			rGraph[u][v] -= path_flow;
			rGraph[v][u] += path_flow;
		}
		max_flow += path_flow;
	}

	// Flow is maximum now, find vertices reachable from s 
	bool visited[V];
	memset(visited, false, sizeof(visited));
	dfs(rGraph, startNode, visited);

	// Print all edges that are from a reachable vertex to 
	// non-reachable vertex in the original graph 
	for (int i = 0; i < V; i++)
		for (int j = 0; j < V; j++)
			if (visited[i] && !visited[j] && graph[i][j])
				cout << i << " - " << j << endl;

	return max_flow;
}

/*
 *	ALGOS
 */
struct uc {
	float ucRatio;
	unsigned int index;
};

bool content[TABSIZE] = {
	true,
	true,
	true,
	true
};

/*
 *	BRANCH & BOUND
 */

struct bbNode {
	unsigned int	lev;
	bool* set;
	float			activity;
	float			cost;
	float			bound;
};

int cmp(const void* elm1, const void* elm2) {
	int a = ((uc*)elm1)->ucRatio;
	int b = ((uc*)elm2)->ucRatio;
	return a > b ? -1 : a < b;
}

int piecewise_linear_distribution_calc(double* cost, double* activity) {
	{
		const int nrolls = 10000; // number of experiments
		const int nstars = 100;   // maximum number of stars to distribute

		std::default_random_engine generator;
		int n = sizeof(TABSIZE) / sizeof(activity[0]);

		std::vector<double> activityVector(activity, activity + TABSIZE);
		boost::random::piecewise_linear_distribution<>
			distribution(&activity[0], &activity[TABSIZE], &cost[0]);

		int p[10] = {};

		for (int i = 0; i < nrolls; ++i) {
			int number = distribution(generator);
			++p[number];
		}

		std::cout << "a piecewise_linear_distribution:" << std::endl;
		for (int i = 0; i < 9; ++i) {
			std::cout << i << "-" << i + 1 << ": ";
			std::cout << std::string(p[i] * nstars / nrolls, '*') << std::endl;
		}
		return 0;
	}
}
float calculateBound(bbNode nd, nodeStruct* node, uc* ucTab, float maxCost, unsigned int size) {
	float addCost = 0.f;
	float addactivity = 0.f;

	for (int i = nd.lev; i < size; i++) {
		if (nd.cost + addCost + node[ucTab[i].index].cost > maxCost) {
			//	Partial add of the best c/u object
			float partialActivity = ((maxCost - addCost) / node[ucTab[i].index].cost * node[ucTab[i].index].duration);
			return nd.activity + addactivity + partialActivity;
		}
		addCost += node[ucTab[i].index].cost;
		addactivity += node[ucTab[i].index].duration;
	}
	return nd.cost + addCost;
}

float calculateActivity(bbNode nd, float* activity, unsigned int size) {
	float res = 0.f;
	for (int i = 0; i < size; i++)
		if (nd.set[i])
			res += activity[i];
	return res;
}

void cpyNode(bbNode* dest, bbNode src, unsigned int size) {
	*dest = src;
	dest->set = (bool*)malloc(size * sizeof(bool));
	memcpy(dest->set, src.set, size * sizeof(bool));
}

void setNextObject(bbNode* nd, bool val, float* cost, uc* ucTab) {
	nd->set[ucTab[nd->lev].index] = val;
	nd->cost += val * cost[ucTab[nd->lev].index];
	nd->lev++;
}

void BranchAndBound(nodeStruct* node, int maxCost, unsigned int size) {
	//	Declarations
	uc* ucTab = (uc*)malloc(size * sizeof(uc));
	queue<bbNode> queueNd;

	//	Sorting items by Cost/activity ratio
	for (int i = 0; i < size; i++) {
		ucTab[i].ucRatio = node[i].duration / node[i].cost;
		ucTab[i].index = i;
	}
	qsort(ucTab, size, sizeof(uc), cmp);

	//	Initializing level, set of items selected etc.
	bbNode bsf;

	bsf.lev = 0;
	bsf.set = (bool*)malloc(size * sizeof(bool));
	bsf.activity = -INFINITY;
	bsf.cost = 0.f;
	bsf.bound = 0.f;
	memset(bsf.set, 0, size * sizeof(bool));

	bbNode aNode;

	aNode.lev = 0;
	aNode.set = (bool*)malloc(size * sizeof(bool));
	aNode.activity = 0.f;
	aNode.cost = 0.f;
	aNode.bound = calculateBound(bsf, node, ucTab, maxCost, size);
	memset(aNode.set, 0, size * sizeof(bool));

	//	Enqueuing
	queueNd.push(aNode);
	int nodeIndex = 0;

	while (!queueNd.empty()) {

		nodeIndex++;
		bbNode currNode = queueNd.front();
		queueNd.pop();

		memcpy(content, currNode.set, size * sizeof(bool));

		float _cost = 0.f;
		for (int i = 0; i < size; i++) {
			_cost += currNode.set[i] * node[i].cost;
		}

		cout << "Cout : " << _cost << endl << endl;

		if (currNode.bound > bsf.activity) {

			//	Create a node nextAdded equal to currNode with the next item added
			bbNode nextAdded;
			cpyNode(&nextAdded, currNode, size);
			delete[] currNode.set;	//	Memory management
			setNextObject(&nextAdded, true, &node[nodeIndex].cost, ucTab);
			nextAdded.bound = calculateBound(nextAdded, &node[nodeIndex], ucTab, maxCost, size);
			nextAdded.activity = calculateActivity(nextAdded, &node[nodeIndex].duration, size);
			//	if nextAdded's activityit > bestSoFar's activityit
			if (nextAdded.activity > bsf.activity) {
				//	set bsf equal to nextAdded
				delete[] bsf.set;
				cpyNode(&bsf, nextAdded, size);
			}

			//	if nextAdded's bound > bestSoFar's activityi
			if (nextAdded.bound > bsf.activity) {
				//	enqueue nextAdded
				queueNd.push(nextAdded);
			}
		}

		//	Create a node nextNotAdded equal to currNode without the next item added
		bbNode nextNotAdded;
		cpyNode(&nextNotAdded, currNode, size);
		delete[] currNode.set;	//	Memory management
		setNextObject(&nextNotAdded, false, &node[nodeIndex].cost, ucTab);

		nextNotAdded.bound = calculateBound(nextNotAdded, &node[nodeIndex], ucTab, maxCost, size);
		nextNotAdded.activity = calculateActivity(nextNotAdded, &node[nodeIndex].duration, size);

		//	if nextNotAdded's bound > bsf's activity
		if (nextNotAdded.bound > bsf.activity) {
			//	enqueue nextNotAdded
			queueNd.push(nextNotAdded);
		}
	}

}

int main() {


	int i, n_activities, top, j;

	std::cout << "############## Critical Path management ################\n\n";
	std::cout << "Enter the number of activities : ";
	cin >> n_activities; // n_activities is the number of activities

	struct nodeStruct* nodes = new nodeStruct[n_activities + 2]; // number of activities here 0th activity is the start
									  // and the (n+1)th activity refers finish both having duration 0

	nodes[0].name = "Start";
	nodes[0].duration = 0;
	nodes[n_activities + 1].name = "Finish";
	nodes[n_activities + 1].duration = 0;
	// input of all the activities
	for (i = 1; i <= n_activities; i++) {
		std::cout << "\n\nEnter activity name #" << i << " : ";
		cin >> nodes[i].name;
		//getline(cin, nodes[i].name);
		std::cout << "Enter duration for " << i << " : ";
		cin >> nodes[i].duration;
		std::cout << "Enter cost for " << i << " : ";
		cin >> nodes[i].duration;
		std::cout << " Enter the es, ef, ls, lf, st" << i << ":";
		cin >> nodes[i].es >> nodes[i].ef >> nodes[i].ls >> nodes[i].lf >> nodes[i].st;

	}
	int startNode, endNode;
	startNode = 0;
	endNode = n_activities - 1;

	std::cout << "\n\n\t\tactivities entered :\n";
	for (i = 0; i <= n_activities + 1; i++) {
		std::cout << "\t\t" << i << ". " << nodes[i].name << " " << nodes[i].duration << endl;
	}


	vector< vector<int> > adj;  // adj represents sucessor list
	vector< vector<int> > pred; // pred reperesents predecessor list
	vector< vector<int> > successorDurations;  // adj represents sucessor list


	// initialization of both lists with empty vectors
	for (i = 0; i <= n_activities; i++) {
		vector<int> temp;
		adj.push_back(temp);
		pred.push_back(temp);
		successorDurations.push_back(temp);
	}

	// initialization of successor list based on user input
	// NOTE : User need to input all the activities with no predecessors as the successor of "Start"
	std::cout << "\n\nNOTE : User need to input all the activities with no predecessors as the successor of \"Start\"";
	for (i = 0; i <= n_activities; i++) {
		std::cout << "\n\nEnter successors for activity " << nodes[i].name << " : ";
		vector<int> temp = ReadNumbers();
		if (temp.size() == 0) {
			adj[i].push_back(n_activities);
			pred[n_activities].push_back(i);
		}
		for (int j = 0; j < temp.size(); j++)
			adj[i].push_back(temp[j]);
		for (int j = 0; j < temp.size(); j++)
			pred[temp[j]].push_back(i);
	}
	//Input Data needs to be given as a set of data
	for (i = 0; i <= n_activities; i++) {
		std::cout << "\n\nEnter successors Duration " << nodes[i].name << " : ";
		vector<int> temp = ReadNumbers();
		for (int j = 0; j < temp.size(); j++)
			successorDurations[i].push_back(temp[j]);
	}

	if (DBG) {

		//debugging
		std::cout << "\nSuccessor matrix :\n";
		for (i = 0; i < n_activities + 1; i++) {
			std::cout << i << " : ";
			for (j = 0; j < adj[i].size(); j++) {
				std::cout << adj[i][j] << "->";
			}
			std::cout << endl;
		}

		std::cout << "Predecessor matrix :\n";
		for (i = 0; i < n_activities + 1; i++) {
			std::cout << i << " : ";
			for (j = 0; j < pred[i].size(); j++) {
				std::cout << pred[i][j] << "->";
			}
			std::cout << endl;
		}
	}

	// calculating earliest start and finish times for each activity
	// topological sort of activity is required here
	stack<int> Stack;
	vector<bool> visit(n_activities + 2, false);
	topologicalSortUtil(0, visit, Stack, adj);

	nodes[0].es = 0;
	nodes[0].ef = 0;
	Stack.pop();

	while (!Stack.empty()) {
		top = Stack.top();
		int max_f = -1;
		for (i = 0; i < pred[top].size(); i++) {
			if (max_f < nodes[pred[top][i]].ef) {
				max_f = nodes[pred[top][i]].ef;
			}
		}
		nodes[top].es = max_f;
		nodes[top].ef = max_f + nodes[top].duration;
		Stack.pop();
	}


	if (DBG) {
		std::cout << "Es and Ef : \n";
		for (i = 0; i < n_activities + 2; i++) {
			std::cout << i << " " << nodes[i].name << " " << nodes[i].es << " " << nodes[i].ef << endl;
		}
	}

	// calculating latest start and finish time for each activity

	stack<int> Stack2;
	vector<bool> visit2(n_activities + 2, false);

	topologicalSortUtil(0, visit2, Stack2, pred);

	nodes[n_activities + 1].ls = nodes[n_activities + 1].es;
	nodes[n_activities + 1].lf = nodes[n_activities + 1].ef;
	Stack2.pop();
	while (!Stack2.empty()) {
		top = Stack2.top();
		int min_s = 99999;
		for (i = 0; i < adj[top].size(); i++) {
			if (min_s > nodes[adj[top][i]].ls) {
				min_s = nodes[adj[top][i]].ls;
			}
		}
		nodes[top].lf = min_s;
		nodes[top].ls = min_s - nodes[top].duration;
		Stack2.pop();
	}

	std::cout << "\n\n";
	if (DBG) {
		std::cout << "Ls and Lf : \n";
		for (i = 0; i < n_activities + 2; i++) {
			std::cout << i << " " << nodes[i].name << " " << nodes[i].ls << " " << nodes[i].lf << endl;
		}
	}


	// display of results
	std::cout << "RESULTS : \n\n";
	std::cout << "\t#\tactivity\tDur.\tEs\tEf\tLs\tLf\tST\n\n";
	for (i = 0; i < n_activities + 2; i++) {
		nodes[i].st = nodes[i].ls - nodes[i].es;
		std::cout << "\t" << i << "\t" << nodes[i].name << "\t" << nodes[i].duration << "\t"
			<< nodes[i].es << "\t" << nodes[i].ef << "\t" << nodes[i].ls << "\t"
			<< nodes[i].lf << "\t" << nodes[i].st << "\n\n";
	}

	//backtracking
	int graph[V][V] = { 0 };
	for (int i = 0; i < V; i++) {
		for (int j = 0; j < V; j++)
			graph[i][j] = successorDurations[i][j];
	}

	// Doing the necessary changes
	int maxDuration = fordFulkerson(graph, n_activities, startNode, endNode);
	std::cout << "The maximum possible flow is " << maxDuration;

	// Branch And Bound algorithm
	BranchAndBound(nodes, maxDuration, n_activities);

	double* durations = new double[n_activities];
	double* cost = new double[n_activities];

	ofstream outputFile("data.txt");

	// initialization of both lists with empty vectors
	for (i = 0; i <= n_activities; i++) {
		durations[i] = nodes[i].duration;
		cost[i] = nodes[i].cost;
		outputFile << i << "\t" << nodes[i].duration << "\t" << std::fixed << nodes[i].cost << "\n";
	}

	outputFile.close();
	plotprocesses();

	//Backtrack and Double check the result
	//PieceWise Linear Distribution
	piecewise_linear_distribution_calc(durations, cost);
	maxDuration = fordFulkerson(graph, n_activities, startNode, endNode);
	std::cout << "The maximum possible flow is " << maxDuration;

	return 0;

}
