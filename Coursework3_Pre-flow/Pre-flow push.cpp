#include<iostream>
#include<stdio.h>
using namespace std;
#define MAX 100

int n;					//vex-number
int H[MAX];				//height
int E[MAX];				//excess
int G[MAX][MAX];		//graph
int R[MAX][MAX];		//residual graph

void push_start(){			//s push-start
	for(int i=0;i<n;i++)
	{
		if(G[0][i]>0)
		{
			E[i]=G[0][i];
			R[0][i]=0;
			R[i][0]=G[0][i];
		}
	}
}
void push(int node){
	for(int i=0;i<n;i++)
	{
		if(R[node][i]>0&&H[node]>H[i]&&E[node]>0)
		{
			int flow=min(R[node][i],E[node]);
			E[i]+=flow;
			E[node]-=flow;
			R[node][i]-=flow;
			R[i][node]+=flow;
		}
	}
}
void relabel(int node){
	int min_H=100;			//bigger number
	for(int i=0;i<n;i++)
	{
		if(R[node][i]>0)
		{
			if(H[i]<min_H)
				min_H=H[i];			
		}
	}
	if(min_H>=H[node])
	{
		H[node]=min_H+1;
	}
}
void preflow_push(){
	int flag=1;
	while(flag)
	{
		flag=0;
		for(int i=1;i<n-1;i++)
		{
			if(E[i]>0)
			{

				flag=1;
				relabel(i);
				push(i);
			}
		}
	}
}
void output(){
	cout<<"capacity:"<<endl;
	for(int i=0;i<n;i++)			
	{
		for(int j=0;j<n;j++)
		{
			if(G[i][j]>0)
			{
				cout<<i+1<<"->";
				cout<<j+1<<" capacity: "<<G[i][j]<<"\t";
			}
		}
		cout<<endl;
	}
	cout<<"maxflow:"<<endl;
	for(int i=0;i<n-1;i++)			
	{
		for(int j=0;j<n;j++)
		{
			if(G[i][j]>0)
			{
				cout<<i+1<<"->";
				cout<<j+1<<" flow: "<<G[i][j]-R[i][j]<<"\t";
			}
		}
		cout<<endl;
	}
}
int main(){
	cin>>n;
	int v1,v2,ca;
	memset(E,0,sizeof(E));
	memset(H,0,sizeof(H));
	memset(G,0,sizeof(G));
	memset(R,0,sizeof(G));
	while(cin>>v1>>v2>>ca&&(v1!=0||v2!=0||ca!=0))
	{
		G[v1][v2]=ca;
		R[v1][v2]=ca;
	}
	H[0]=n;
	push_start();
	preflow_push();
	output();
	system("pause");
	return 0;
}