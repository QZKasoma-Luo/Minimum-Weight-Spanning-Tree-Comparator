/* PrimVsKruskal.java
   CSC 226 - Fall 2023
   Assignment 3 - Prim MST versus Kruskal MST Template
   
   The file includes the "import edu.princeton.cs.algs4.*;" so that you can use
   any of the code in the algs4.jar file. You should be able to compile your program
   with the command
   
	javac -cp .;algs4.jar PrimVsKruskal.java
	
   To conveniently test the algorithm with a large input, create a text file
   containing a test graphs (in the format described below) and run
   the program with
   
	java -cp .;algs4.jar PrimVsKruskal file.txt
	
   where file.txt is replaced by the name of the text file. Note: different operating systems have different commands.
   You should all know how to run the algs4.jar file from lab work.
   
   The input consists of a graph (as an adjacency matrix) in the following format:
   
    <number of vertices>
	<adjacency matrix row 1>
	...
	<adjacency matrix row n>
	
   Entry G[i][j] >= 0.0 of the adjacency matrix gives the weight (as type double) of the edge from 
   vertex i to vertex j (if G[i][j] is 0.0, then the edge does not exist).
   Note that since the graph is undirected, it is assumed that G[i][j]
   is always equal to G[j][i].
*/

 import edu.princeton.cs.algs4.*;

import java.util.HashSet;
import java.util.Scanner;
 import java.io.File;

//Do not change the name of the PrimVsKruskal class
public class PrimVsKruskal{

	//Add and inital the variables from algs4.jar, some of the variables and comments are from alg4.jar
	private EdgeWeightedGraph EW_Graph; //transform the input txt file to a edge weighted graph
	private Queue<Edge> mst = new Queue<>();
	private IndexMinPQ<Double> pq; //eligible crossing edges
	private boolean[] marked; //true if v in MST
	private Edge[] edgeTo;        // edgeTo[v] = shortest edge from tree vertex to non-tree vertex
    private double[] distTo;      // distTo[v] = weight of shortest such edge

	//Transfor the input txt file to a edge weighted graph variable to be able to pass to the eagerPrim and Kruskal function
	public void parameterTransfrom(double[][] G){
		this.EW_Graph = new EdgeWeightedGraph(G.length);
		for(int i = 0; i < G.length; i++){
			for(int j = 0; j < G.length; j++){
				if(G[i][j] > 0.0){
					Edge e = new Edge(i, j, G[i][j]);
					EW_Graph.addEdge(e);
				}
			}
		}
		marked = new boolean[EW_Graph.V()];
		pq = new IndexMinPQ<>(EW_Graph.V());
	}
	
	public void eagerPrim(EdgeWeightedGraph G){
		//initialize the edgeTo and distTo
		edgeTo = new Edge[G.V()];
		distTo = new double[G.V()];
		for(int v = 0; v < G.V(); v++){
			distTo[v] = Double.POSITIVE_INFINITY;
		}
		
		//initialize the first vertex
		distTo[0] = 0.0; //Assume G is connected
		pq.insert(0,0.0);
		while(!pq.isEmpty()){
			visit(G, pq.delMin()); //repeatedly delete the min weight edge e = v-w from pq 
		}
	}

	private void visit(EdgeWeightedGraph G, int v){
		marked[v] = true;
		for(Edge e:G.adj(v)){ // add v to tree; update data structures
			int w = e.other(v); // for each edge e = v-w, add w to pq if not already in tree
			if(marked[w]) continue;
			if(e.weight()<distTo[w]){
				edgeTo[w] = e; //add e to tree
				distTo[w] = e.weight();
				if(pq.contains(w)) pq.changeKey(w, distTo[w]); // update distance to w or insert distance to w
				else pq.insert(w, distTo[w]);
			}
		}
	}

	public Iterable<Edge> edges(){ //create the MST
		mst = new Queue<Edge>();
		for (int v = 0; v < edgeTo.length; v++){
			Edge e = edgeTo[v];
			if(e != null){
				mst.enqueue(e);
			}
		}
		return mst;
	}	
	/* PrimVsKruskal(G)
		Given an adjacency matrix for connected graph G, with no self-loops or parallel edges,
		determine if the minimum spanning tree of G found by Prim's algorithm is equal to 
		the minimum spanning tree of G found by Kruskal's algorithm.
		
		If G[i][j] == 0.0, there is no edge between vertex i and vertex j
		If G[i][j] > 0.0, there is an edge between vertices i and j, and the
		value of G[i][j] gives the weight of the edge.
		No entries of G will be negative.
	*/

	static boolean PrimVsKruskal(double[][] G){
		int n = G.length;

		/* Build the MST by Prim's and the MST by Kruskal's */
		/* (You may add extra methods if necessary) */
		
		/* ... Your code here ... */
		
		PrimVsKruskal pmTransfer = new PrimVsKruskal();
		pmTransfer.parameterTransfrom(G);
		pmTransfer.eagerPrim(pmTransfer.EW_Graph);
		Iterable<Edge> pmMST = pmTransfer.edges();


		EdgeWeightedGraph kruskalGraph = new EdgeWeightedGraph(n);
		for(int i = 0; i < n; i++){
			for(int j = i+1; j < n; j++){
				if(G[i][j] > 0.0){
					Edge e = new Edge(i, j, G[i][j]);
					kruskalGraph.addEdge(e);
				}
			}
		}

		KruskalMST kruskal_MST = new KruskalMST(kruskalGraph);
		Iterable<Edge> kruskalMST = kruskal_MST.edges();
		/* Determine if the MST by Prim equals the MST by Kruskal */
		boolean pvk = true;
		/* ... Your code here ... */

		HashSet<Edge> pmSet = new HashSet<>();
		for(Edge e: pmMST){
			pmSet.add(e);
		}

		HashSet<Edge> kruskalSet = new HashSet<>();
		for(Edge e: kruskalMST){
			kruskalSet.add(e);
		}

		if(!pmSet.equals(kruskalSet)){
			pvk = false;
		}

		return pvk;	
	}
	
		
	/* main()
	   Contains code to test the PrimVsKruskal function. You may modify the
	   testing code if needed, but nothing in this function will be considered
	   during marking, and the testing process used for marking will not
	   execute any of the code below. 
	*/
   public static void main(String[] args) {
		Scanner s;
		if (args.length > 0){
			try{
				s = new Scanner(new File(args[0]));
			} catch(java.io.FileNotFoundException e){
				System.out.printf("Unable to open %s\n",args[0]);
				return;
			}
			System.out.printf("Reading input values from %s.\n",args[0]);
		}else{
			s = new Scanner(System.in);
			System.out.printf("Reading input values from stdin.\n");
		}
		
		int n = s.nextInt();
		double[][] G = new double[n][n];
		int valuesRead = 0;
		for (int i = 0; i < n && s.hasNextDouble(); i++){
			for (int j = 0; j < n && s.hasNextDouble(); j++){
				G[i][j] = s.nextDouble();
				if (i == j && G[i][j] != 0.0) {
					System.out.printf("Adjacency matrix contains self-loops.\n");
					return;
				}
				if (G[i][j] < 0.0) {
					System.out.printf("Adjacency matrix contains negative values.\n");
					return;
				}
				if (j < i && G[i][j] != G[j][i]) {
					System.out.printf("Adjacency matrix is not symmetric.\n");
					return;
				}
				valuesRead++;
			}
		}
		
		if (valuesRead < n*n){
			System.out.printf("Adjacency matrix for the graph contains too few values.\n");
			return;
		}	
		
        boolean pvk = PrimVsKruskal(G);
        System.out.printf("Does Prim MST = Kruskal MST? %b\n", pvk);
    }
}
