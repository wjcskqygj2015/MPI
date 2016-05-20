#include <mpi/mpi.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <ctime>
#include <sys/time.h>
#include <immintrin.h>
#include <stdlib.h>
#include <vector>
const int maxN=1<<10;
float A[maxN][maxN];
float B[maxN][maxN];
float C[maxN][maxN];
double error=0.0;
int totalProcNum=0;
int nowProcID=0;
MPI::Datatype array;
void init(int n)
{
	srand(time(NULL));
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			C[i][j] = B[i][j]= A[i][j] = (float)(rand()%100+1);
}

void LUResolveSeq(int n,float resolving[][maxN])
{
	for(int k=0;k<n;k++)
    {
        for(int j=k+1;j<n;j++)
            resolving[k][j]/=resolving[k][k];
        //resolving[k][k]=1;
        for(int i=k+1;i<n;i++)
            for(int j=k+1;j<n;j++)
                resolving[i][j]=resolving[i][j]-(resolving[i][k]*resolving[k][j]);
        //    resolving[i][k]=0;
    }

}

float buffer[maxN]={};
void LUResolveMPI(int n,float resolving[][maxN])
{
//First We need to Matrix resolving, for it is common to initailize resolving in One process
	array=MPI::FLOAT.Create_contiguous(n);
	array.Commit();

//Then we need to Scatter the content of resolving
	for(int i=0;i<n;i+=totalProcNum)
		MPI::COMM_WORLD.Scatter(resolving+i, 1, array, resolving+i+nowProcID, 1, array, 0);

//This time we just assume that ProcNum can completely divide the Scale, then we finished 

//We start the round	

	for(int nowRound=0;nowRound<n;nowRound++)
	{
		int master=nowRound%totalProcNum;
		if(nowProcID==master)
			for(int i=nowRound+1;i<n;i++)
				buffer[i]=(resolving[nowRound][i]/=resolving[nowRound][nowRound]);
		//Then we use broadcast let this part to all	
		MPI::COMM_WORLD.Bcast(buffer+nowRound+1,n-nowRound-1,MPI::INT,master);
		//Now we can use the buffer Number to be calculate the number
		for(int row=n-totalProcNum+nowProcID;row>=nowRound+1;row-=totalProcNum)
			for(int column=nowRound+1;column<n;column++)
				resolving[row][column]=resolving[row][column]-resolving[row][nowRound]*buffer[column];
	}
//We finished calculate process

//Then we gather the whole matrix
	for(int i=0;i<n;i+=totalProcNum)
		MPI::COMM_WORLD.Gather(resolving+i+nowProcID,1,array,resolving+i,1,array,0);
//We finish gathering matrix
	array.Free();
}


void LUResolveMPIPipeLine(int n,float resolving[][maxN])
{

}


void detection()
{
    error=0;
    for(int i=0;i<maxN;i++)
        for(int j=i+1;j<maxN;j++)
            error+=fabs(A[i][j]-B[i][j]);

    //for(int i=0;i<maxN;i++)
    //    for(int j=i+1;j<maxN;j++)
    //        error+=fabs(A[i][j]-C[i][j]);
}


int main(int argc, char **argv)
{
	long long L1,L2;
	timeval tv;
//Initialize MPI enviroment
	MPI::Init(argc,argv);
	totalProcNum=MPI::COMM_WORLD.Get_size();
	nowProcID=MPI::COMM_WORLD.Get_rank();	
	
	if(nowProcID==0)
	{
		init(maxN);
		gettimeofday(&tv,NULL);
		L1=tv.tv_sec*1000*1000+tv.tv_usec;
		LUResolveSeq(maxN,A);
		gettimeofday(&tv,NULL);
		L2=tv.tv_sec*1000*1000+tv.tv_usec;
		printf("Sequen time: %lldus\n",L2-L1);	
	}
	MPI::COMM_WORLD.Barrier();
//There just use barrier to sychronize and prevent other Processes begin to timing


	gettimeofday(&tv,NULL);
	L1=tv.tv_sec*1000*1000+tv.tv_usec;
	LUResolveMPI(maxN,B);
	gettimeofday(&tv,NULL);
	L2=tv.tv_sec*1000*1000+tv.tv_usec;
	printf("I am Process %d, MPINop time: %lldus\n",nowProcID,L2-L1);		

//	gettimeofday(&tv,NULL);
//	L1=tv.tv_sec*1000*1000+tv.tv_usec;
	//LUResolveMPIPipeLine(maxN,C);
//	gettimeofday(&tv,NULL);
//	L2=tv.tv_sec*1000*1000+tv.tv_usec;
//	printf("I am Process %d, MPIPip time: %lldus\n",nowProcID,L2-L1);		
	
	if(nowProcID==0)
	{
		detection();
		printf("errors is: %lf\n",error);
	}	
	MPI::Finalize();
	return 0;

}



