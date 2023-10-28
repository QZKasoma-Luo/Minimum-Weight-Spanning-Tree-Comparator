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
import java.util.Scanner;
import java.io.File;


//Do not change the name of the PrimVsKruskal class
public class PrimVsKruskal {

	/*
	 * PrimVsKruskal(G)
	 * Given an adjacency matrix for connected graph G, with no self-loops or
	 * parallel edges,
	 * determine if the minimum spanning tree of G found by Prim's algorithm is
	 * equal to
	 * the minimum spanning tree of G found by Kruskal's algorithm.
	 * 
	 * If G[i][j] == 0.0, there is no edge between vertex i and vertex j
	 * If G[i][j] > 0.0, there is an edge between vertices i and j, and the
	 * value of G[i][j] gives the weight of the edge.
	 * No entries of G will be negative.
	 */

	 enum EdgeList{
		inMst,
		excluded,
		unknown
	 }

	static boolean PrimVsKruskal(double[][] G){
		int n = G.length;

		/* Build the MST by Prim's and the MST by Kruskal's */
		/* (You may add extra methods if necessary) */

		/* ... Your code here ... */

		IndexMinPQ<Edge> primPQ = new IndexMinPQ<Edge>(n); // choose the lightest edge and delete it from the pq and add it to mst
		double[] primDistTo = new double[n]; //the shortest edge from the current vertices to the mst
		int primEdgeTo[] = new int[n]; // store the parent vertices of the current vertices
		boolean[] primMarked = new boolean[n]; // ture if the vertices is in the Prim mst
		Queue<Edge> primMst = new Queue<Edge>();
		EdgeList[][] edgeStatus = new EdgeList[n][n];


		IndexMinPQ<Edge> kruskalPQ = new IndexMinPQ<Edge>(n * ((n) / 2)); //the krusukal will load all the edge to the pq at the beginning, that we need to set the possible max size of the graph
		Queue<Edge> kruskalMst = new Queue<Edge>();
        Edge[] kruskalEdges = new Edge[n * ((n-1) / 2)]; //the krusukal will load all the edge to the pq at the beginning, that we need to set the possible max size of the graph
		UF kruskalUF = new UF(n); //Union find to balance the mst

		// intialize the prim
		for (int veretx = 0; veretx < n; veretx++) { 
				primDistTo[veretx] = Double.POSITIVE_INFINITY; //set all the edge weight to infinity to indicate the mst edge connection status.
		}
		primDistTo[0] = 0.0;
		primPQ.insert(0, new Edge(0, 0, 0.0)); // add the starting vertex 0 with no any edges to the primPQ
		for(int i = 0; i < n; i++){
			for(int k = 0; k< n; k++){
				edgeStatus[i][k] = EdgeList.unknown;
			}
		}

		// intial the kruskal
		int kruskalIndex = 0; //use to store the edges
		for(int i = 0; i < n; i++){
			for(int j = i + 1; j < n; j++){
				if(G[i][j] > 0.0){ //check if veretx i and veretx j has an edge that connects each other
					Edge edge = new Edge(i, j, G[i][j]); // create an Edge object from the G
					kruskalEdges[kruskalIndex] = edge; //insert all the edge to the kruskalEdgeList
					kruskalPQ.insert(kruskalIndex, edge); // insert the Edge object into the kruskalPQ
					kruskalIndex++;
				}
			}
		}

		/* Determine if the MST by Prim equals the MST by Kruskal */
		boolean pvk = true;
		/* ... Your code here ... */
		while(!primPQ.isEmpty() || !kruskalPQ.isEmpty()){
			
			//eagerPrim part
			if(!primPQ.isEmpty()){
				int minEdgeVertex = primPQ.delMin(); //deletes and return the veretx associated with the smallest edge
				primMarked[minEdgeVertex] = true; //add the vertext to the mst
				if(minEdgeVertex != 0){ //we initialed the primPQ with starting veretx 0, so we need to make sure we don't repeatly add the vertex 0 again
					Edge edge = new Edge(primEdgeTo[minEdgeVertex], minEdgeVertex, G[primEdgeTo[minEdgeVertex]][minEdgeVertex]);
					primMst.enqueue(edge);
				}

				for(int p = 0; p<n; p++){ 
					if(G[minEdgeVertex][p] > 0 && !primMarked[p]){
						if(G[minEdgeVertex][p] < primDistTo[p]){
							primDistTo[p] = G[minEdgeVertex][p];
							primEdgeTo[p] = minEdgeVertex;
							if(primPQ.contains(p)){
								primPQ.changeKey(p, new Edge(p, p, primDistTo[p]));
							}else{
								primPQ.insert(p, new Edge(p, p, primDistTo[p]));
							}
							edgeStatus[minEdgeVertex][p] = EdgeList.inMst; // Marking as included
							edgeStatus[p][minEdgeVertex] = EdgeList.inMst; // Marking as included
							
						} else {
								edgeStatus[minEdgeVertex][p] = EdgeList.excluded; 
								edgeStatus[p][minEdgeVertex] = EdgeList.excluded; 
						}
					}
				}
			}	
			
			//kruskal part
			if(!kruskalPQ.isEmpty()){
				int edgeIndex = kruskalPQ.delMin();
				Edge e = kruskalEdges[edgeIndex];
				int v = e.either();
				int m = e.other(v);
				if(kruskalUF.find(v) != kruskalUF.find(m)){
					kruskalUF.union(v, m);
					kruskalMst.enqueue(e);
				}
				if(edgeStatus [v][m] == EdgeList.excluded){
					pvk = false;
				}
			}
		}	
		double primTotalWeight = 0.0;
		System.out.println("Prim MST:");
		for (Edge e : primMst) {
			int v = e.either();
			int w = e.other(v);
			System.out.println(v + " - " + w + " : " + e.weight());
			primTotalWeight += e.weight();
		}
		System.out.println("Total weight of Prim MST: " + primTotalWeight);
		double kruskalTotalWeight = 0.0;
		
		System.out.println("Kruskal MST:");
		for (Edge e : kruskalMst) {
			int v = e.either();
			int w = e.other(v);
			System.out.println(v + " - " + w + " : " + e.weight());
			kruskalTotalWeight += e.weight();
		}
			System.out.println("Total weight of Kruskal MST: " + kruskalTotalWeight);

			return pvk;	
	}

	/*
	 * main()
	 * Contains code to test the PrimVsKruskal function. You may modify the
	 * testing code if needed, but nothing in this function will be considered
	 * during marking, and the testing process used for marking will not
	 * execute any of the code below.
	 */
	public static void main(String[] args) {
		Scanner s;
		if (args.length > 0) {
			try {
				s = new Scanner(new File(args[0]));
			} catch (java.io.FileNotFoundException e) {
				System.out.printf("Unable to open %s\n", args[0]);
				return;
			}
			System.out.printf("Reading input values from %s.\n", args[0]);
		} else {
			s = new Scanner(System.in);
			System.out.printf("Reading input values from stdin.\n");
		}

		int n = s.nextInt();
		double[][] G = new double[n][n];
		int valuesRead = 0;
		for (int i = 0; i < n && s.hasNextDouble(); i++) {
			for (int j = 0; j < n && s.hasNextDouble(); j++) {
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

		if (valuesRead < n * n) {
			System.out.printf("Adjacency matrix for the graph contains too few values.\n");
			return;
		}

		boolean pvk = PrimVsKruskal(G);
		System.out.printf("Does Prim MST = Kruskal MST? %b\n", pvk);
	}
}
