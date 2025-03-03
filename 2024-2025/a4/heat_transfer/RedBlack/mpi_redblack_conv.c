#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include "utils.h"
#define TEST_CONV
typedef enum { false, true } bool;


void Jacobi(double ** u_previous, double ** u_current, int X_min, int X_max, int Y_min, int Y_max) {
        int i,j;
        for (i=X_min;i<X_max;i++)
                for (j=Y_min;j<Y_max;j++)
                        u_current[i][j]=(u_previous[i-1][j]+u_previous[i+1][j]+u_previous[i][j-1]+u_previous[i][j+1])/4.0;
}

void GaussSeidel(double ** u_previous, double ** u_current, int X_min, int X_max, int Y_min, int Y_max, double omega) {
        int i,j;
        for (i=X_min;i<X_max;i++)
                for (j=Y_min;j<Y_max;j++)
                        u_current[i][j]=u_previous[i][j]+(u_current[i-1][j]+u_previous[i+1][j]+u_current[i][j-1]+u_previous[i][j+1]-4*u_previous[i][j])*omega/4.0;
}

void RedSOR(double ** u_previous, double ** u_current, int X_min, int X_max, int Y_min, int Y_max, double omega) {
        int i,j;
        for (i=X_min;i<X_max;i++)
                for (j=Y_min;j<Y_max;j++)
                        if ((i+j)%2==0)
                                u_current[i][j]=u_previous[i][j]+(omega/4.0)*(u_previous[i-1][j]+u_previous[i+1][j]+u_previous[i][j-1]+u_previous[i][j+1]-4*u_previous[i][j]);
}

void BlackSOR(double ** u_previous, double ** u_current, int X_min, int X_max, int Y_min, int Y_max, double omega) {
        int i,j;
        for (i=X_min;i<X_max;i++)
                for (j=Y_min;j<Y_max;j++)
                        if ((i+j)%2==1)
                                u_current[i][j]=u_previous[i][j]+(omega/4.0)*(u_current[i-1][j]+u_current[i+1][j]+u_current[i][j-1]+u_current[i][j+1]-4*u_previous[i][j]);
}


int main(int argc, char ** argv) {
    int rank,size;
    int global[2],local[2]; //global matrix dimensions and local matrix dimensions (2D-domain, 2D-subdomain)
    int global_padded[2];   //padded global matrix dimensions (if padding is not needed, global_padded=global)
    int grid[2];            //processor grid dimensions
    int i,j,t;
    int global_converged=0,converged=0; //flags for convergence, global and per process
    MPI_Datatype dummy;     //dummy datatype used to align user-defined datatypes in memory
    double omega; 			//relaxation factor - useless for Jacobi

    struct timeval tts,ttf,tcs,tcf;   //Timers: total-> tts,ttf, computation -> tcs,tcf
    double ttotal=0,tcomp=0,total_time,comp_time;
    
    double ** U, ** u_current, ** u_previous, ** swap; //Global matrix, local current and previous matrices, pointer to swap between current and previous
    

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    //----Read 2D-domain dimensions and process grid dimensions from stdin----//

    if (argc!=5) {
        fprintf(stderr,"Usage: mpirun .... ./exec X Y Px Py");
        exit(-1);
    }
    else {
        global[0]=atoi(argv[1]);
        global[1]=atoi(argv[2]);
        grid[0]=atoi(argv[3]);
        grid[1]=atoi(argv[4]);
    }

    //----Create 2D-cartesian communicator----//
	//----Usage of the cartesian communicator is optional----//

    MPI_Comm CART_COMM;         //CART_COMM: the new 2D-cartesian communicator
    int periods[2]={0,0};       //periods={0,0}: the 2D-grid is non-periodic
    int rank_grid[2];           //rank_grid: the position of each process on the new communicator
		
    MPI_Cart_create(MPI_COMM_WORLD,2,grid,periods,0,&CART_COMM);    //communicator creation
    MPI_Cart_coords(CART_COMM,rank,2,rank_grid);	                //rank mapping on the new communicator

    //----Compute local 2D-subdomain dimensions----//
    //----Test if the 2D-domain can be equally distributed to all processes----//
    //----If not, pad 2D-domain----//
    
    for (i=0;i<2;i++) {
        if (global[i]%grid[i]==0) {
            local[i]=global[i]/grid[i];
            global_padded[i]=global[i];
        }
        else {
            local[i]=(global[i]/grid[i])+1;
            global_padded[i]=local[i]*grid[i];
        }
    }

	//Initialization of omega
    omega=2.0/(1+sin(3.14/global[0]));

    //----Allocate global 2D-domain and initialize boundary values----//
    //----Rank 0 holds the global 2D-domain----//
    if (1) {
        U=allocate2d(global_padded[0],global_padded[1]);   
        init2d(U,global[0],global[1]);
    }

    //----Allocate local 2D-subdomains u_current, u_previous----//
    //----Add a row/column on each size for ghost cells----//

    u_previous=allocate2d(local[0]+2,local[1]+2);
    u_current=allocate2d(local[0]+2,local[1]+2);   
       
    //----Distribute global 2D-domain from rank 0 to all processes----//
         
 	//----Appropriate datatypes are defined here----//
	/*****The usage of datatypes is optional*****/
    
    //----Datatype definition for the 2D-subdomain on the global matrix----//

    MPI_Datatype global_block;
    MPI_Type_vector(local[0],local[1],global_padded[1],MPI_DOUBLE,&dummy);
    MPI_Type_create_resized(dummy,0,sizeof(double),&global_block);
    MPI_Type_commit(&global_block);

    //----Datatype definition for the 2D-subdomain on the local matrix----//

    MPI_Datatype local_block;
    MPI_Type_vector(local[0],local[1],local[1]+2,MPI_DOUBLE,&dummy);
    MPI_Type_create_resized(dummy,0,sizeof(double),&local_block);
    MPI_Type_commit(&local_block);

    //----Rank 0 defines positions and counts of local blocks (2D-subdomains) on global matrix----//
    int * scatteroffset, * scattercounts;
    if(1) {
        scatteroffset=(int*)malloc(size*sizeof(int));
        scattercounts=(int*)malloc(size*sizeof(int));
        for (i=0;i<grid[0];i++)
            for (j=0;j<grid[1];j++) {
                scattercounts[i*grid[1]+j]=1;
                scatteroffset[i*grid[1]+j]=(local[0]*local[1]*grid[1]*i+local[1]*j);
            }
    }


    //----Rank 0 scatters the global matrix----//
    
    //----Rank 0 scatters the global matrix----//

	//*************TODO*******************//



	/*Fill your code here*/

	

	/*Make sure u_current and u_previous are
		both initialized*/

        zero2d(u_current,local[0]+2,local[1]+2);
        zero2d(u_previous,local[0]+2,local[1]+2);
        
	MPI_Scatterv(*U,scattercounts,scatteroffset,global_block,(*u_current+local[1]+3),1,local_block,0,CART_COMM);	
        for(i = 0; i < local[0]+2; i++){
                for(j = 0; j < local[1]+2; j++){
                        u_previous[i][j]=u_current[i][j];
                }
        }



     //************************************//


    if (rank==0)
        free2d(U);

 
     
	//----Define datatypes or allocate buffers for message passing----//

	//*************TODO*******************//



	/*Fill your code here*/

    MPI_Datatype row;
    MPI_Type_contiguous(local[1],MPI_DOUBLE,&dummy);
    MPI_Type_create_resized(dummy,0,sizeof(double),&row);
    MPI_Type_commit(&row);

    MPI_Datatype column;
    MPI_Type_vector(local[0],1,local[1]+2,MPI_DOUBLE,&dummy);
    MPI_Type_create_resized(dummy,0,sizeof(double),&column);
    MPI_Type_commit(&column);








	//************************************//


    //----Find the 4 neighbors with which a process exchanges messages----//

	//*************TODO*******************//
    int north, south, east, west;


	/*Fill your code here*/


	/*Make sure you handle non-existing
		neighbors appropriately*/
   // printf("%d ",rank);
    MPI_Cart_shift(CART_COMM, 0, 1, &north, &south);
    MPI_Cart_shift(CART_COMM, 1,1, &west, &east);
/*
    printf("I am rank %d and I found a north neighbor : %d\n",rank,north);
    printf("I am rank %d and I found a south neighbor : %d\n",rank,south);
    printf("I am rank %d and I found a west neighbor : %d\n",rank,west);
    printf("I am rank %d and I found a east neighbor : %d\n",rank,east);
*/
   // printf("I am rank %d and my coords are (%d,%d)\n",rank,rank_grid[0],rank_grid[1]);

	//************************************//



    //---Define the iteration ranges per process-----//
	//*************TODO*******************//

    int i_min,i_max,j_min,j_max;
    int padding[2];
	
    padding[1]=global_padded[1]-global[1];
    padding[0]=global_padded[0]-global[0];
	/*Fill your code here*/


	bool is_internal(int x, int y){
		return x!=0 && y!=0 && x!=(grid[0]-1) && y!=(grid[1]-1);
	}

	bool is_right_padded(){
	     return global[1]%grid[1]!=0;
	}
	
        bool is_bottom_padded(){
             return global[0]%grid[0]!=0;
        }

	
	bool is_boundary_not_padded(int x, int y){
		return !is_internal(x,y) && !is_right_padded() && !is_bottom_padded();
	}

	if(is_internal(rank_grid[0],rank_grid[1])){
		i_min=1; // skip ghost cells
		j_min=1;
		i_max=local[0]+1; //max value -> i<i_max not equal in for loop	
		j_max=local[1]+1;
	}
	else if(is_boundary_not_padded(rank_grid[0],rank_grid[1])){
		if(north<0 && west<0 && south>=0 && east>=0){ //north-west border
			i_min=2;
			j_min=2;
			i_max=local[0]+1;
			j_max=local[1]+1;
		}
		else if(west<0 && north>=0 && east>=0 && south>=0){ //middle west borders
			i_min=1;
			j_min=2;
			i_max=local[0]+1;
			j_max=local[1]+1;
		}
		else if(west<0 && south<0 && east>=0 && north>=0){ // south west border
			i_min=1;
			j_min=2;
			i_max=local[0];
			j_max=local[1]+1;
		}
		else if(south<0 && west>=0 && north>=0 && east>=0){ // south middle borders
			i_min=1;
			j_min=1;
			i_max=local[0];
			j_max=local[1]+1;
		}
		else if(east<0 && south<0 && north>=0 && west>=0){ // south east border
			i_min=1;
			j_min=1;
			i_max=local[0];
			j_max=local[1];
		}
		else if(east<0 && north>=0 && west>=0 && south>=0){ //east middle borders
			i_min=1;
			j_min=1;
			i_max=local[0]+1;
			j_max=local[1];
		}
		else if(north<0 && east<0 && west>=0 && south>=0){ //north east border
			i_min=2;
			j_min=1;
			i_max=local[0]+1;
			j_max=local[1];
		}
		else if(north<0 && east>=0 && south>=0 && west>=0){ //middle north borders
			i_min=2;
			j_min=1;
			i_max=local[0]+1;
			j_max=local[1]+1;
		}
		else if(west<0 && north<0 && south<0 && east>=0){//we have only east neighbor
			i_min=2;
			j_min=2;
			i_max=local[0];
			j_max=local[1]+1;
		}
		else if(east<0 && north<0 && south<0 && west>=0){ //only west neighbor
			i_min=2;
			j_min=1;
			i_max=local[0];
			j_max=local[1];
		}
		else if(east<0 && north<0 && south<0 && west<0){ //no neighbors (I am a lonely sunflower)
			i_min=2;
			j_min=2;
			i_max=local[0];
			j_max=local[1];
		}
	}
	
	else if(is_right_padded() && is_bottom_padded()){
                if(north<0 && east<0 && south>=0 && west>=0){ //north-east border
                        i_min=2;
                        j_min=1;
                        i_max=local[0]+1;
                        j_max=local[1]-padding[1];
                }
                else if(east<0 && north>=0 && south>=0 && west>=0){ //middle east borders
                        i_min=1;
                        j_min=1;
                        i_max=local[0]+1;
                        j_max=local[1]-padding[1];
                }
                else if(east<0 && south<0 && north>=0 && west>=0){ //south east border
                        i_min=1;
                        j_min=1;
                        i_max=local[0]-padding[0];
                        j_max=local[1]-padding[1];
                }
                else if(south<0 && north>=0 && east>=0 && west>=0){ //middle south borders
                        i_min=1;
                        j_min=1;
                        i_max=local[0]-padding[0];
                        j_max=local[1]+1;
                }
                else if(south<0 && west<0 && north>=0 && east>=0){ //south west border
                        i_min=1;
                        j_min=2;
                        i_max=local[0]-padding[0];
                        j_max=local[1]+1;
                }

	}
	else if(is_right_padded()){
                if(north<0 && east<0 && south>=0 && west>=0){ //north-east border
                        i_min=2;
                        j_min=1;
                        i_max=local[0]+1;
                        j_max=local[1]-padding[1];
                }
                else if(east<0 && north>=0 && south>=0 && west>=0){ //middle east borders
                        i_min=1;
                        j_min=1;
                        i_max=local[0]+1;
                        j_max=local[1]-padding[1];
                }
                else if(east<0 && south<0 && north>=0 && west>=0){ //south east border
                        i_min=1;
                        j_min=1;
                        i_max=local[0];
                        j_max=local[1]-padding[1];
                }
	        else if(west<0 && north<0 && south<0 && east>=0){//we have only east neighbor
                        i_min=2;
                        j_min=2;
                        i_max=local[0];
                        j_max=local[1]+1;
                }
                else if(east<0 && north<0 && south<0 && west>=0){ //only west neighbor
                        i_min=2;
                        j_min=1;
                        i_max=local[0];
                        j_max=local[1]-padding[1];
                }

		
	}
	
	else if(is_bottom_padded()){
                if(east<0 && south<0 && north>=0 && west>=0){ //south east border
                        i_min=1;
                        j_min=1;
                        i_max=local[0]-padding[0];
                        j_max=local[1];
                }
                else if(south<0 && north>=0 && east>=0 && west>=0){ //middle south borders
                        i_min=1;
                        j_min=1;
                        i_max=local[0]-padding[0];
                        j_max=local[1]+1;
                }
                else if(south<0 && west<0 && north>=0 && east>=0){ //south west border
                        i_min=1;
                        j_min=2;
                        i_max=local[0]-padding[0];
                        j_max=local[1]+1;
                }

	}


	
	//printf("I am rank %d with coords=(%d, %d) and my boundaries are: i_min=%d, i_max=%d, j_min=%d, j_max=%d\n",rank,rank_grid[0],rank_grid[1],i_min,i_max,j_min,j_max);
	/*Three types of ranges:
		-internal processes
		-boundary processes
		-boundary processes and padded global array
	*/





	//************************************//



	MPI_Barrier(MPI_COMM_WORLD);
 	//----Computational core----//   
	gettimeofday(&tts, NULL);
    #ifdef TEST_CONV
    for (t=0;t<T && !global_converged;t++) {
    #endif
    #ifndef TEST_CONV
    #undef T
    #define T 256
    for (t=0;t<T;t++) {
    #endif


	 	//*************TODO*******************//
     
 

//	printf("Got into computational core, t=%d, rank=%d\n",t,rank);


		/*Fill your code here*/


		/*Compute and Communicate*/

		/*Add appropriate timers for computation*/
	MPI_Request requests[8];
	int counter = 0;	

	if(north>=0){
		//MPI_Sendrecv(*u_current+local[1]+3,1,row,north,0,*u_current+1,1,row,north,0,MPI_COMM_WORLD,&requests[counter]);

		MPI_Isend(*u_current+local[1]+3,1,row,north,0,MPI_COMM_WORLD,&requests[counter]);
		counter+=1;
		MPI_Irecv(*u_current+1,1,row,north,0,MPI_COMM_WORLD,&requests[counter]);
		counter+=1;
	}
	if(south>=0){
		//MPI_Sendrecv(*u_current+local[0]*(local[1]+2)+1,1,row,south,1,*u_current+((local[1]+2)*(local[0]+1))+1,1,row,south,1,MPI_COMM_WORLD,&requests[counter]);
	
		MPI_Isend(*u_current+local[0]*(local[1]+2)+1,1,row,south,0,MPI_COMM_WORLD,&requests[counter]);
		counter+=1;
		MPI_Irecv(*u_current+((local[1]+2)*(local[0]+1))+1,1,row,south,0,MPI_COMM_WORLD,&requests[counter]);
		counter+=1;
	}
	if(east>=0){
		//MPI_Sendrecv(*u_current+2*(local[1]+1),1,column,east,3,*u_current+2*(local[1]+1)+1,1,column,east,3,MPI_COMM_WORLD,&requests[counter]);
		
		MPI_Isend(*u_current+2*(local[1]+1),1,column,east,0,MPI_COMM_WORLD,&requests[counter]);
		counter+=1;
		MPI_Irecv(*u_current+2*(local[1]+1)+1,1,column,east,0,MPI_COMM_WORLD,&requests[counter]);
		counter+=1;
	}
	if(west>=0){
		//MPI_Sendrecv(*u_current+local[1]+3,1,column,west,4,*u_current+local[1]+2,1,column,west,4,MPI_COMM_WORLD,&requests[counter]);
	
		MPI_Isend(*u_current+local[1]+3,1,column,west,0,MPI_COMM_WORLD,&requests[counter]);
		counter+=1;
		MPI_Irecv(*u_current+local[1]+2,1,column,west,0,MPI_COMM_WORLD,&requests[counter]);
		counter+=1;
	}
	
//	printf("Before waitall, rank=%d, t=%d\n",rank,t);
	MPI_Status status[8];
	MPI_Waitall(counter,requests,status);

	
	swap=u_previous;
	u_previous=u_current;
	u_current=swap;	

	gettimeofday(&tcs,NULL);	
	RedSOR(u_previous,u_current,i_min,i_max,j_min,j_max,omega);
	gettimeofday(&tcf,NULL);
	tcomp+=(tcf.tv_sec-tcs.tv_sec)+(tcf.tv_usec-tcs.tv_usec)*0.000001;
	MPI_Request Brequests[8];
	int Bcounter = 0;	

	if(north>=0){
		MPI_Isend(*u_current+local[1]+3,1,row,north,0,MPI_COMM_WORLD,&Brequests[Bcounter]);
		Bcounter+=1;
		MPI_Irecv(*u_current+1,1,row,north,0,MPI_COMM_WORLD,&Brequests[Bcounter]);
		Bcounter+=1;
	}
	if(south>=0){
		MPI_Isend(*u_current+local[0]*(local[1]+2)+1,1,row,south,0,MPI_COMM_WORLD,&Brequests[Bcounter]);
		Bcounter+=1;
		MPI_Irecv(*u_current+((local[1]+2)*(local[0]+1))+1,1,row,south,0,MPI_COMM_WORLD,&Brequests[Bcounter]);
		Bcounter+=1;
	}
	if(east>=0){
		MPI_Isend(*u_current+2*(local[1]+1),1,column,east,0,MPI_COMM_WORLD,&Brequests[Bcounter]);
		Bcounter+=1;
		MPI_Irecv(*u_current+2*(local[1]+1)+1,1,column,east,0,MPI_COMM_WORLD,&Brequests[Bcounter]);
		Bcounter+=1;
	}
	if(west>=0){
		MPI_Isend(*u_current+local[1]+3,1,column,west,0,MPI_COMM_WORLD,&Brequests[Bcounter]);
		Bcounter+=1;
		MPI_Irecv(*u_current+local[1]+2,1,column,west,0,MPI_COMM_WORLD,&Brequests[Bcounter]);
		Bcounter+=1;
	}
	
//	printf("Before waitall, rank=%d, t=%d\n",rank,t);
	MPI_Waitall(Bcounter,Brequests,status);


	gettimeofday(&tcs,NULL);
	//printf("Welcome to Jacobi Mr rank %d, t=%d\n",rank,t);
	//Jacobi(u_previous,u_current,i_min,i_max,j_min,j_max);
	//GaussSeidel(u_previous, u_current,i_min, i_max, j_min, j_max, omega);
	BlackSOR(u_previous,u_current,i_min,i_max,j_min,j_max,omega);
	gettimeofday(&tcf,NULL);
	tcomp+=(tcf.tv_sec-tcs.tv_sec)+(tcf.tv_usec-tcs.tv_usec)*0.000001;
//	printf("Bye bye Jacobi from rank %d, t=%d",rank,t);




		#ifdef TEST_CONV
        if (t%C==0) {
			//*************TODO**************//
  			/*Test convergence*/
		converged=converge(u_previous,u_current,i_min,i_max-1,j_min,j_max-1);		
	
		MPI_Allreduce(&converged, &global_converged, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
		
		}		
		#endif


		//************************************//        
    }

    MPI_Barrier(MPI_COMM_WORLD);
    gettimeofday(&ttf,NULL);

    ttotal=(ttf.tv_sec-tts.tv_sec)+(ttf.tv_usec-tts.tv_usec)*0.000001;

    MPI_Reduce(&ttotal,&total_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&tcomp,&comp_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);



    //----Rank 0 gathers local matrices back to the global matrix----//
   
    if (rank==0) {
            U=allocate2d(global_padded[0],global_padded[1]);
    }


	//*************TODO*******************//


	
	/*Fill your code here*/

	MPI_Gatherv(*u_current+local[1]+3,1,local_block,*U,scattercounts,scatteroffset,global_block,0,CART_COMM);



     //************************************//


	
   

	//----Printing results----//

	//**************TODO: Change "Jacobi" to "GaussSeidelSOR" or "RedBlackSOR" for appropriate printing****************//
    if (rank==0) {
        printf("GaussSeidelSOR X %d Y %d Px %d Py %d Iter %d ComputationTime %lf TotalTime %lf midpoint %lf\n",global[0],global[1],grid[0],grid[1],t,comp_time,total_time,U[global[0]/2][global[1]/2]);
	
        #ifdef PRINT_RESULTS
        char * s=malloc(50*sizeof(char));
        sprintf(s,"resJacobiMPI_%dx%d_%dx%d",global[0],global[1],grid[0],grid[1]);
        fprint2d(s,U,global[0],global[1]);
        free(s);
        #endif

    }
    MPI_Finalize();
    return 0;
}
