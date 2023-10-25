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
public class PrimVsKruskal{

	/* PrimVsKruskal(G)
		Given an adjacency matrix for connected graph G, with no self-loops or parallel edges,
		determine if the minimum spanning tree of G found by Prim's algorithm is equal to 
		the minimum spanning tree of G found by Kruskal's algorithm.
		
		If G[i][j] == 0.0, there is no edge between vertex i and vertex j
		If G[i][j] > 0.0, there is an edge between vertices i and j, and the
		value of G[i][j] gives the weight of the edge.
		No entries of G will be negative.
	*/
	private static Queue<double[]> eagerPrim(double[][] G, int G_length){
		double[] distTo = new double[G_length]; //the shortest edge from the current vertices to the mst
		int edgeTo[] = new int[G_length]; // store the parent vertices of the current vertices
		boolean[] marked = new boolean[G_length]; // ture if the vertices is in the mst
		IndexMinPQ<Double> pq = new IndexMinPQ<Double>(G_length);// choose the lightest edge and delete it from the pq and add it to mst
		Queue<double[]> primMst = new Queue<double[]>(); // store the mst

		for(int i = 0; i < G_length; i++){
			distTo[i] = Double.POSITIVE_INFINITY;
		}
		distTo[0] = 0.0;
		pq.insert(0, 0.0);

		while(!pq.isEmpty()){ // while pq is not empty it means there are still vertices not in the mst
			int minEdgeVertex = pq.delMin(); // delete the lightest edge to a veretx from pq and return that veretx to minEdgeVertex
			marked[minEdgeVertex] = true; // mark the current minEdgeVertex is visited
			if(minEdgeVertex != 0){ // check if the vertex is not the root
				// add the edge to the mst
				primMst.enqueue(new double[] {edgeTo[minEdgeVertex], minEdgeVertex, G[edgeTo[minEdgeVertex]][minEdgeVertex]});
			}

			for(int p = 0; p < G_length; p++){ // check all the vertices to update their shortest edge to the mst
				if(G[minEdgeVertex][p] > 0.0 && !marked[p]){ //if the edge from minEdgeVertex to p is not 0 and p is not in the mst
					if(G[minEdgeVertex][p] < distTo[p]){ // if the edge from minEdgeVertex to p is lighter than the current shortest edge to the mst
						distTo[p] = G[minEdgeVertex][p]; // update the shortest edge to the mst
						edgeTo[p] = minEdgeVertex; // update the parent vertices of p
						if(pq.contains(p)){ // if p is in the pq, update the shortest edge to the mst
							pq.changeKey(p, distTo[p]);
						}else{ // if p is not in the pq, add it to the pq
							pq.insert(p, distTo[p]);
						}
					}
				}
			}
		}
		return primMst;
	}

	private static Queue<double[]> kruskal(double[][] G, int G_length){
		Queue<double[]> krusualMst = new Queue<double[]>();
		Edge[] edges = new Edge[G_length * ((G_length-1) / 2)];
		int index = 0;
		
		for(int i = 0; i < G_length; i++){
			for(int j = i+1; j < G_length; j++){
				if(G[i][j] > 0.0){
					edges[index++] = new Edge(i, j, G[i][j]);
				}
			}
		}

		Array.sort(edges, Comparator.comparingDouble(Edge::weight));
	}

	static boolean PrimVsKruskal(double[][] G){
		int n = G.length;

		/* Build the MST by Prim's and the MST by Kruskal's */
		/* (You may add extra methods if necessary) */

		/* ... Your code here ... */
		Queue<double[]> primResult = eagerPrim(G, n);
    	System.out.println("Prim Tree:");
		double totalWeight = 0.0;
		while (!primResult.isEmpty()) {
			double[] edge = primResult.dequeue();
			System.out.printf("%.0f-%.0f %.5f\n", edge[0], edge[1], edge[2]);
			totalWeight += edge[2];
		}
		
		System.out.printf("%.5f\n", totalWeight);
		/* Determine if the MST by Prim equals the MST by Kruskal */
		boolean pvk = true;
		/* ... Your code here ... */

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
