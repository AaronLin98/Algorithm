#include<iostream>
#include<stdio.h>
#include<string.h>
using namespace std;
#define MAXINT 32767
#define MVNUM 100
typedef struct{
	char vexs[MVNUM];
	int arcs[MVNUM][MVNUM];
	int vexnum,arcsnum;
}AMGraph;
void output(AMGraph G,int w,int path[]){
	if(path[w]==-1) 
	{
		cout<<G.vexs[w];
	}
	else
	{
		int e=path[w];
		output(G,e,path);
		cout<<" ";
		cout<<G.vexs[w];
	}
}
int locatevex(AMGraph G,char v1){
	int i=0;
	for(i=0;i<G.vexnum;i++)
	{	
		if(v1==G.vexs[i])
			return i;
	}
}
void input(AMGraph &G,int n,int m){
	int i=0,j=0,k=0;
	char v1,v2;
	int w=0;
	G.vexnum=n;
	G.arcsnum=m;
	for(i=0;i<G.vexnum;i++)
		cin>>G.vexs[i];
	for(i=0;i<G.vexnum;i++)
		for(j=0;j<G.vexnum;j++)
			{	
				G.arcs[i][j]=MAXINT;
			}
	for(k=0;k<G.arcsnum;k++)
	{
		cin>>v1>>v2>>w;
		i=locatevex(G,v1);
		j=locatevex(G,v2);
		G.arcs[i][j]=w;
	}
}
void shortestpath_DIJ(AMGraph G,int v0,int ve){
	int S[MVNUM];
	int D[MVNUM];
	int path[MVNUM];
	int min=0;
	int n=G.vexnum;
	int v=0,i=0,w=0;
	for(v=0;v<n;v++)
	{
		S[v]=false;
		D[v]=G.arcs[v0][v];
		if(D[v]<MAXINT) path[v]=v0;
		else path[v]=-1;
	}
	S[v0]=true;
	D[v0]=0;					//初始化
	for(i=1;i<n;i++)
	{
		min=MAXINT;
		for(w=0;w<n;w++)
			if(!S[w]&&D[w]<min)
			{v=w;min=D[w];}
		S[v]=true;				//选择最短路径
		if(S[ve]==true) 
		{
			cout<<D[ve]<<endl;
			output(G,ve,path);
			cout<<endl;
			break;
		}
		for(w=0;w<n;w++)
			if(!S[w]&&(D[v]+G.arcs[v][w]<D[w]))
			{	
				D[w]=D[v]+G.arcs[v][w];
				path[w]=v;
			}					//更新最短路径
	}

}
int main(){
	int n=0,m=0;
	AMGraph G;
	char a,b;
	int v0=0,ve=0;
	cin>>n>>m;
	while(n!=0&&m!=0){
		input(G,n,m);
		cin>>a>>b;
		v0=locatevex(G,a);
		ve=locatevex(G,b);
		shortestpath_DIJ(G,v0,ve);
		cin>>n>>m;
	}
	return 0;
}
