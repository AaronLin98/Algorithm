//151303120 林哲宇
#include<iostream>
#include<stdio.h>
#include<limits.h>
using namespace std;
#define MAX_INT 1999	//无穷大值
#define MAX 20

void Graph(int G[][MAX],int n){
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			G[i][j]=MAX_INT;
		}
	}
	cout<<"请输入各条边的大小,当输入为0 0 0时结束"<<endl;
	int a,b,w;
	while(cin>>a&&cin>>b&&cin>>w)
	{
		if(a==0&&b==0&&w==0)
		{
			break;
		}
		else
		{
			G[a][b]=w;
		}
	}
}
void M_initialition(int M[][MAX],int n){
	M[0][n-1]=0;
	for(int i=0;i<n-1;i++)	
	{
		M[0][i]=MAX_INT;
	}
}
int min_int(int a,int b){
	if(a>b)
		return b;
	else
		return a;
}
int min_M(int G[][MAX],int M[][MAX],int v,int k,int n){
	int w;
	int min;
	for(int i=0;i<n;i++)
	{
		if(i==0)
		{
			min=G[v][i]+M[k][i];
			w=i;
		}
		else
		{
			if(G[v][i]+M[k][i]<min)
			{
				min=G[v][i]+M[k][i];
				w=i;
			}
		}
	}
	return w;
}
void Shortest_Path(int G[][MAX],int M[][MAX],int n){
	for(int i=1;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			int w=min_M(G,M,j,i-1,n);	//寻找最小边所在的点
			M[i][j]=min_int(M[i-1][j],M[i-1][w]+G[j][w]);
		}
	}
}
void output(int M[][MAX],int k,int n){
	for(int i=0;i<n;i++)
	{
		cout<<M[k-1][i]<<endl;
	}
}
void out_put(int M[MAX][MAX],int n){
	cout<<"对于OPT(n-1,V)而言结果为："<<endl;
	output(M,n-1,n);
	cout<<"对于OPT(n,V)而言结果为："<<endl;
	output(M,n,n);
}
int judge(int M[MAX][MAX],int n){
	for(int i=0;i<n;i++)
	{
		if(M[n-2][i]!=M[n-1][i])
		{
			return 0;
		}
	}
	return 1;
}
void out_result(int M[MAX][MAX],int n){
	int flag=judge(M,n);
	if(flag)
	{
		cout<<"******图中无负环******\n"<<endl;
	}
	else 
	{
		cout<<"******图中有负环******\n"<<endl;
	}
}

int main(){
	int n;
	cout<<"请输入图中点的个数"<<endl;
	cin>>n;
	int M[MAX][MAX];			//初始化
	M_initialition(M,n);
	int G[MAX][MAX];
	Graph(G,n);
	Shortest_Path(G,M,n);		//计算
	out_put(M,n);				//输出
	out_result(M,n);
	system("pause");
	return 0;
}