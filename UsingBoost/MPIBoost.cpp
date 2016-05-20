#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <ctime>
#include <sys/time.h>
#include <immintrin.h>
#include <stdlib.h>
#include <mpi/mpi.h>
using namespace std;

const int maxN=1<<10;

float A[maxN][maxN];
float B[maxN][maxN];
float C[maxN][maxN];
float X[maxN][maxN];
double error=0.0;

//进行全局层面的初始化，由于使用的是构造函数，所以可以保证是在所有的其他函数执行之前进行初始化

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


int sendbuffer[maxN]={};
int recvbuffer[maxN]={};
void LUResolveMPI(int n,float resolving[][maxN])
{
	int totalProcNum=0;
	int nowProcId=0;
	 
	for(int nowRound=0;nowRound<n;nowRound++)
	{
		int master=nowRound%totalProcNum;
		if(nowProcId==master)
		{
			for(int i=nowRound+1;i<n;i++)
				sendbuffer[i]=(resolving[nowRound][i]/=resolving[nowRound][nowRound]);
		}
	}		
/*
	if (world.rank() == 0) 
	{
		mpi::request reqs[2];
		gps_position msg, out_msg = gps_position(1,2,3);
    reqs[0] = world.isend(1, 0, out_msg);
    reqs[1] = world.irecv(1, 1, msg);
    mpi::wait_all(reqs, reqs + 2);
    std::cout << msg << "!" << std::endl;
  } 
	else 
	{
    mpi::request reqs[2];
    gps_position msg, out_msg = gps_position(3,2,1);
    reqs[0] = world.isend(0, 1, out_msg);
    reqs[1] = world.irecv(0, 0, msg);
    mpi::wait_all(reqs, reqs + 2);
    std::cout << msg << ", ";
  }
*/
}

void LUResolveMPIPipeLine(int n,float resolving[][maxN])
{

/*
	if (world.rank() == 0) 
	{
		mpi::request reqs[2];
		gps_position msg, out_msg = gps_position(1,2,3);
    reqs[0] = world.isend(1, 0, out_msg);
    reqs[1] = world.irecv(1, 1, msg);
    mpi::wait_all(reqs, reqs + 2);
    std::cout << msg << "!" << std::endl;
  } 
	else 
	{
    mpi::request reqs[2];
    gps_position msg, out_msg = gps_position(3,2,1);
    reqs[0] = world.isend(0, 1, out_msg);
    reqs[1] = world.irecv(0, 0, msg);
    mpi::wait_all(reqs, reqs + 2);
    std::cout << msg << ", ";
  }
*/

}

void detection()
{
    error=0;
    for(int i=0;i<maxN;i++)
        for(int j=i+1;j<maxN;j++)
            error+=fabs(A[i][j]-B[i][j]);

    for(int i=0;i<maxN;i++)
        for(int j=i+1;j<maxN;j++)
            error+=fabs(A[i][j]-C[i][j]);
}



int main()
{


	long long L1,L2;
	timeval tv;
	init(maxN);
	gettimeofday(&tv,NULL);
	L1=tv.tv_sec*1000*1000+tv.tv_usec;
	LUResolveSeq(maxN,A);
	gettimeofday(&tv,NULL);
	L2=tv.tv_sec*1000*1000+tv.tv_usec;
	printf("Sequen time: %lldus\n",L2-L1);	
	
	gettimeofday(&tv,NULL);
	L1=tv.tv_sec*1000*1000+tv.tv_usec;
	LUResolveMPI(maxN,B);
	gettimeofday(&tv,NULL);
	L2=tv.tv_sec*1000*1000+tv.tv_usec;
	printf("MPINop time: %lldus\n",L2-L1);		

	gettimeofday(&tv,NULL);
	L1=tv.tv_sec*1000*1000+tv.tv_usec;
	//LUResolveMPIPipeLine(maxN,C);
	gettimeofday(&tv,NULL);
	L2=tv.tv_sec*1000*1000+tv.tv_usec;
	printf("MPIPip time: %lldus\n",L2-L1);

	detection();
	printf("absoluteErr: %lf\n",error);

	return 0;
}
