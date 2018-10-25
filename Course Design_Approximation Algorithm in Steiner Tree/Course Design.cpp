#include<iostream>
#include<cmath>
#include<vector>
#include<set>
#include<queue>
#include<fstream>
#include<algorithm>
using namespace std;

#define MAX 100
#define PI 3.1416

typedef struct POINT //顶点 
{
	double x;
	double y;
	int index;

	friend bool operator < (struct POINT const &v1, struct POINT const &v2)
	{
		if (v1.x == v2.x)
		{
			return v1.y < v2.y;
		}
		else
		{
			return v1.x < v2.x;
		}
	}

	friend bool operator == (struct POINT const &v1, struct POINT const &v2)
	{
		return (v1.x == v2.x&&v1.y == v2.y);
	}

	friend bool operator !=(const POINT &a, const POINT &b)  
	{
		if (a.x == b.x && a.y == b.y)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

}V;

typedef struct EDGE {
	V a;
	V b;
	int index;
	double cost;
	vector<int> MST_T_index;	//所在的MST三角形的序列号 0<=MST_T_index.size()<=2

	friend bool operator < (struct EDGE const &e1, struct EDGE const &e2)
	{
		if (e1.a == e2.a)
		{
			return e1.b < e2.b;
		}
		else
		{
			return e1.a < e2.a;
		}
	}

	friend bool operator == (struct EDGE const &e1, struct EDGE const &e2)
	{
		return (e1.a == e2.a&&e1.b == e2.b) || (e1.a == e2.b&&e1.b == e2.a);
	}

}E;

typedef struct TRIANGLE 
{
	POINT a;
	POINT b;
	POINT c;
	int index;
	int MST_index;	//在MST里面的序列号 

	friend bool operator ==(const TRIANGLE &t1, const TRIANGLE &t2) //判断两个三角形是否相等 
	{
		if (t1.a != t2.a &&t1.a != t2.b &&t1.a != t2.c)
			return false;
		else if (t1.b != t2.a &&t1.b != t2.b &&t1.b != t2.c)
			return false;
		else if (t1.c != t2.a &&t1.c != t2.b &&t1.c != t2.c)
			return false;
		else
			return true;
	}

}T;

struct ACS {		//相邻关系
	int index;							//MST三角形的本身的序列号index
	int ab_index, bc_index, ac_index;	//MST三角形的ab,bc,ac边所相邻的MST三角形的序列号(若无则为-1) 
};					

struct edge			//记录相邻三角形的边
{
	POINT a, b;
	int adjtriangle;
};

struct stPOINT		//斯坦纳点
{
	int index;
	POINT point;
};

struct triasmt		//计算和比较smt需要用
{
	double smt;
	int index;
	POINT stainer;
};


//三角剖分
int n; //点的总数 
vector<POINT> points; //平面上的点集 
vector<TRIANGLE> triangles; //已确定的三角形列表
vector<TRIANGLE> temp_triangles; //未确定的三角形列表 
vector<TRIANGLE> split_triangles; //待分割的三角形列表 
//MST
double G[MAX][MAX];
set<V> MST_set;
set<E> edgesets;
vector<E> edges;
vector<E> MST_edges;
//adjecent
vector<ACS> acs;
int PT_edges[MAX][MAX];
vector<T> MST_triangles;
//生成斯坦纳点
vector<ACS> exisarc(acs);
vector<POINT> steiner;
vector<stPOINT> stainer;
//SMT / MST
double mst_cost = 0;
double smt_cost = 0;


//三角剖分
bool Compare(POINT a, POINT b)	//按横坐标对点进行排序比较
{
	if (a.x <= b.x)
	{
		if (a.x == b.x && a.y > b.y)
			return false;
		else
			return true;
	}
	else
		return false;
}
bool CompareY(POINT a, POINT b)	//按纵坐标对点进行排序比较
{
	if (a.y <= b.y)
	{
		if (a.y == b.y && a.x > b.x)
			return false;
		else
			return true;
	}
	else
		return false;
}
void Input() //输入平面点集及其坐标并排序 
{
	POINT temp;

	cout << "请输入平面点集中顶点数：";
	cin >> n;
	for (int i = 0; i < n; i++)
	{
		cout << "请输入第" << i + 1 << "个点的x、y坐标：";
		cin >> temp.x >> temp.y;
		points.push_back(temp);
	}

	sort(points.begin(), points.end(), Compare);
	//index
	for (int i = 0; i<points.size(); i++)
	{
		points[i].index = i;
	}
}
TRIANGLE Makesupert() //确定超级三角形 
{
	TRIANGLE super;
	double left, right, up, down;	//分别找到点集在四个方向最边缘的坐标值
	double height, length;	//标出该四点组成的方形的长和高
	POINT a, b, c;

	left = points.front().x;
	right = points.front().x;
	up = points.front().y;
	down = points.front().y;
	for (vector<POINT>::iterator it = points.begin(); it != points.end(); it++)
	{
		if (left > it->x)
			left = it->x;
		if (right < it->x)
			right = it->x;
		if (up < it->y)
			up = it->y;
		if (down > it->y)
			down = it->y;
	}

	height = up - down;
	length = right - left;

	a.x = left + length / 2;	//此时求出的三个坐标组成的三角形恰好包围住点集
	a.y = up + length / 2;
	b.x = left - height;
	b.y = down;
	c.x = right + height;
	c.y = down;

	a.y += height / 2;
	b.x -= height;
	b.y -= height / 2;
	c.x += height;
	c.y -= height / 2;

	super.a = a;
	super.b = b;
	super.c = c;

	return super;
}
double Distance(POINT a, POINT b) //求两点间距离 
{
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}
double GetR(TRIANGLE t) //求出三角形外接圆半径 
{
	double R;
	double e1, e2, e3;
	e1 = Distance(t.a, t.b);
	e2 = Distance(t.b, t.c);
	e3 = Distance(t.c, t.a);
	R = e1*e2*e3 / sqrt((e1 + e2 + e3)*(e2 + e3 - e1)*(e1 - e2 + e3)*(e1 + e2 - e3));

	return R;
}
POINT GetPoint(TRIANGLE t) //求出三角形外接圆圆心坐标 
{
	POINT temp;
	temp.x = (t.a.x + t.b.x + t.c.x) / 3;
	temp.y = (t.a.y + t.b.y + t.c.y) / 3;
	return temp;
}
int PosCompare(TRIANGLE t, POINT p) //比较点相对于三角形外接圆的位置 
{
	double R = GetR(t);
	POINT o = GetPoint(t);
	if (Distance(p, o) > R)
	{
		if (p.x > o.x + R)
			return 1; //点在外接圆右侧，为delaunay三角形
		else
			return 0; //点在外接圆外侧，待确定 
	}
	else
		return -1; //点在外接圆内侧，非delaunay三角形 

}
bool Existed(TRIANGLE t)	//判断该三角形是否已存在于不确定列表中
{
	for (vector<TRIANGLE>::iterator it = temp_triangles.begin(); it != temp_triangles.end(); it++)
	{
		if (t == *it)
			return true;
		else
			continue;
	}
	return false;

}
void Split(POINT p)	//将非delaunay三角形分割并去重存入temp_triangles
{
	TRIANGLE temp;
	for (vector<TRIANGLE>::iterator it = split_triangles.begin(); it != split_triangles.end(); it++)
	{
		temp.a = it->a;
		temp.b = it->b;
		temp.c = p;
		if (!Existed(temp))
			temp_triangles.push_back(temp);
		temp.a = it->b;
		temp.b = it->c;
		temp.c = p;
		if (!Existed(temp))
			temp_triangles.push_back(temp);
		temp.a = it->c;
		temp.b = it->a;
		temp.c = p;
		if (!Existed(temp))
			temp_triangles.push_back(temp);
	}
}
bool Related(TRIANGLE super, TRIANGLE t) //检查当前三角形是否与超级三角形有关 
{
	if (t.a != super.a && t.a != super.b && t.a != super.c)
	{
		if (t.b != super.a && t.b != super.b && t.b != super.c)
		{
			if (t.c != super.a && t.c != super.b && t.c != super.c)
				return false;
			else
				return true;
		}
		else
			return true;
	}
	else
		return true;
}
void Delaunay() //Delaunay三角剖分 
{
	int judge;
	int i = 1;
	TRIANGLE super = Makesupert();
	temp_triangles.push_back(super); //将超级三角形存入未确定列表

	for (vector<POINT>::iterator it = points.begin(); it != points.end(); it++) //按顺序遍历所有点 
	{
		for (vector<TRIANGLE>::iterator iter = temp_triangles.begin(); iter != temp_triangles.end();) //遍历未确定列表中的三角形 
		{
			judge = PosCompare(*iter, *it); //判断三角形是否delaunay三角形或不确定
			if (judge == 1)
			{
				triangles.push_back(*iter); //是delaunay三角形，存入确定列表 
				iter = temp_triangles.erase(iter);
			}
			else if (judge == 0)
			{
				iter++; //无法确定，跳过 
			}

			else if (judge == -1)
			{
				split_triangles.push_back(*iter);
				iter = temp_triangles.erase(iter);
			}
		}
		Split(*it);
	}
	//点的遍历结束，合并triangles与temp_triangles 
	triangles.insert(triangles.end(), temp_triangles.begin(), temp_triangles.end());
	//除去与超级三角形有关的三角形 
	for (vector<TRIANGLE>::iterator it = triangles.begin(); it != triangles.end(); )
	{
		if (Related(super, *it))
			it = triangles.erase(it);
		else
			it++;
	}
	//此时得到的triangles列表中为所有三角剖分中的三角形 
	cout << triangles.size() << endl;
	//index
	for (int i = 0; i<triangles.size(); i++)
	{
		triangles[i].index = i;
	}
	//output
	for (vector<TRIANGLE>::iterator it = triangles.begin(); it != triangles.end(); it++)
	{
		cout << "(" << it->a.x << "," << it->a.y << ")  (" << it->b.x << "," << it->b.y << ")  (" << it->c.x << "," << it->c.y << ")" << endl;
	}
}
//MST
double Compute_cost(V v1, V v2) {
	double cost = sqrt((v1.x - v2.x)*(v1.x - v2.x) + (v1.y - v2.y)*(v1.y - v2.y));
	return cost;
}
bool cmp(E e1, E e2) {
	return e1.cost<e2.cost;
}
E reverse(E e) {
	E temp;
	temp.a = e.b;
	temp.b = e.a;
	temp.cost = e.cost;
	return temp;
}
void join_edges(T t) {
	E e1, e2, e3;
	double c1, c2, c3;
	c1 = Compute_cost(t.a, t.b);
	c2 = Compute_cost(t.c, t.b);
	c3 = Compute_cost(t.a, t.c);
	e1.a = t.a;
	e1.a.index = t.a.index;
	e1.b = t.b;
	e1.b.index = t.b.index;
	e1.cost = c1;
	if (!edgesets.count(e1) && !edgesets.count(reverse(e1)))
	{
		edgesets.insert(e1);
		edges.push_back(e1);
	}
	e2.a = t.c;
	e2.a.index = t.c.index;
	e2.b = t.b;
	e2.b.index = t.b.index;
	e2.cost = c2;
	if (!edgesets.count(e2) && !edgesets.count(reverse(e2)))
	{
		edgesets.insert(e2);
		edges.push_back(e2);
	}
	e3.a = t.c;
	e3.a.index = t.c.index;
	e3.b = t.a;
	e3.b.index = t.a.index;
	e3.cost = c3;
	if (!edgesets.count(e3) && !edgesets.count(reverse(e3)))
	{
		edgesets.insert(e3);
		edges.push_back(e3);
	}
}
void MST() {
	memset(G, -1, sizeof(G));
	for (vector<T>::iterator it = triangles.begin(); it != triangles.end(); it++)
	{
		T t1 = (*it);
		join_edges(t1);
	}
	sort(edges.begin(), edges.end(), cmp);
	for (int i = 0; i<edges.size(); i++)
	{
		edges[i].index = i;
	}
	for (vector<E>::iterator it = edges.begin(); it != edges.end(); it++)
	{
		E e = (*it);
		if (!MST_set.count(e.a) || !MST_set.count(e.b))
		{
			MST_edges.push_back(e);
			MST_set.insert(e.a);
			MST_set.insert(e.b);
			G[e.a.index][e.b.index] = e.cost;
			G[e.b.index][e.a.index] = e.cost;
			mst_cost+=e.cost;
		}
	}
/*	for (int i = 0; i<n; i++)
	{
		for (int j = 0; j<n; j++)
		{
			if (G[i][j]>0)
			{
				mst_cost += G[i][j];
			}
		}
	}
	mst_cost = mst_cost / 2.0;*/
	cout << "\nmst:\n" << "mst cost : " << mst_cost  << endl;

}
//寻找相邻三角形
int judge_mst_triangle(V a, V b) {			//判断mst里面两条相邻的边是否能形成一个三角形
	E e;
	e.a = a;
	e.b = b;
	vector<E>::iterator result = find(edges.begin(), edges.end(), e);
	if (result != edges.end())
		return 1;
	else
		return 0;
}
int Triangles_index(V a, V b, V c) {	//得到原三角形的序列号
	T t;
	t.a = a;
	t.b = b;
	t.c = c;
	vector<T>::iterator Tresult = find(triangles.begin(), triangles.end(), t);
	if (Tresult == triangles.end()) 
		return -1;
	int index_T = (*Tresult).index;
	return index_T;
}										
										//push操作
void MST_edge_triangles_push(E e, int index_T) {
	vector<E>::iterator it = find(edges.begin(), edges.end(), e);
	(*it).MST_T_index.push_back(index_T);
}
										//添加MST中的边所属的三角形
void MST_edge_triangles(int index_T) {
	E ab, bc, ac;
	ab.a = triangles[index_T].a;
	ab.b = triangles[index_T].b;

	MST_edge_triangles_push(ab, index_T);


	bc.a = triangles[index_T].b;
	bc.b = triangles[index_T].c;

	MST_edge_triangles_push(bc, index_T);

	ac.a = triangles[index_T].a;
	ac.b = triangles[index_T].c;

	MST_edge_triangles_push(ac, index_T);
}
int ACS_index(E e, int index_T) {				//ACS辅助操作：得到一条边所邻的三角形
	vector<E>::iterator qt = find(edges.begin(), edges.end(), e);
	auto p = find(MST_edges.begin(), MST_edges.end(), (*qt));
	if (p == MST_edges.end()) return -2;
	else if ((*qt).MST_T_index.size() == 2)
	{
		if ((*qt).MST_T_index[0] == index_T)
			return (*qt).MST_T_index[1];
		else
			return (*qt).MST_T_index[0];
	}
	else return -1;
		
}												//ACS辅助操作：得到原三角形的边
void edges_for_triangle(E &ab, E &bc, E &ac, int index_T) {
	ab.a = triangles[index_T].a;
	ab.b = triangles[index_T].b;

	ac.a = triangles[index_T].a;
	ac.b = triangles[index_T].c;

	bc.a = triangles[index_T].b;
	bc.b = triangles[index_T].c;
}
void adjcent() {
	int MST_index = 0;			//记录MST三角形在MST三角形向量中的序列号
	for (int i = 0; i<n; i++)
	{
		int k = 0;
		V a, b, c;
		for (int j = 0; j<n; j++)
		{
			if (G[i][j]>0)
			{
				if (!k)
				{
					k++;
					a = points[j];
				}
				else
				{
					k = 0;
					b = points[j];
					if (judge_mst_triangle(a, b))	//判断是否在原三角剖分中存在该三角形
					{
						c = points[i];				
						int index_T = Triangles_index(a, b, c);			//寻找对应的原来的三角形的序号
						if (index_T == -1) continue;				
						triangles[index_T].MST_index = MST_index;		//形成MST三角形的序列号
						MST_triangles.push_back(triangles[index_T]);	//所有可能形成的MST三角形
						MST_index++;
						MST_edge_triangles(index_T);					//存储每条边所有可能相邻的三角形 （最多两个）
					}
				}
			}
		}
	}
	for (vector<T>::iterator it = MST_triangles.begin(); it != MST_triangles.end(); it++)
	{
		E ab, bc, ac;
		struct ACS acs1;
		int index_T = (*it).index;

		acs1.index = index_T;

		edges_for_triangle(ab, bc, ac, index_T);
		acs1.ab_index = ACS_index(ab, index_T);
		acs1.bc_index = ACS_index(bc, index_T);
		acs1.ac_index = ACS_index(ac, index_T);

		acs.push_back(acs1);
	}
	cout<<"\nMST三角形:"<<endl;
	for (vector<ACS>::iterator it = acs.begin(); it != acs.end(); it++)
	{
		cout << " " << (*it).index << " " << (*it).ab_index << " " << (*it).ac_index << " " << (*it).bc_index << endl;
	}
	cout<<endl;
}
//生成斯坦纳点
double fabsarea(POINT p1, POINT p2, POINT p3)//求面积
{
	return fabs((p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x));
}
POINT getcross(POINT a1, POINT a2, POINT b1, POINT b2)//求线段交点 公式 
{
	double s1, s2;
	s1 = fabsarea(a1, a2, b1);
	s2 = fabsarea(a1, a2, b2);
	POINT t;
	t.x = (b2.x*s1 + b1.x*s2) / (s1 + s2); t.y = (b2.y*s1 + b1.y*s2) / (s1 + s2);
	return t;
}
double stainerpoint(TRIANGLE t, POINT&stainer, int i)
{
	POINT q, w, e;//w max
	EDGE ab, bc, ca;
	ab.cost = pow(t.a.x - t.b.x, 2) + pow(t.a.y - t.b.y, 2); 
	ab.cost = sqrt(ab.cost); 
	ab.a = t.a;
	ab.b = t.b;
	bc.cost = pow(t.b.x - t.c.x, 2) + pow(t.b.y - t.c.y, 2);
	bc.cost = sqrt(bc.cost); bc.a = t.b;
	bc.b = t.c;
	ca.cost = pow(t.c.x - t.a.x, 2) + pow(t.c.y - t.a.y, 2); 
	ca.cost = sqrt(ca.cost); 
	ca.a = t.c;
	ca.b = t.a;
	auto p = find(MST_edges.begin(), MST_edges.end(), ab); //找到三角形在mst中的两条边 
	if (p == MST_edges.end())
	{
		w = t.c; 
		q = t.a;
		e = t.b;
	}
	 p = find(MST_edges.begin(), MST_edges.end(), bc); 
	 if (p == MST_edges.end()) 
	 {
		w = t.a;
		q = t.b; 
		e = t.c;
	}
	 p = find(MST_edges.begin(), MST_edges.end(), ca);
	 if (p == MST_edges.end()) 
	 {
		w = t.b; 
		q = t.c;
		e = t.a;
	}

	POINT trq, trq1, trq2, tre, tre1;
	trq.x = q.x - w.x; trq.y = q.y - w.y;
	//构造等边三角形 
	double j1, j2;
	j1 = 1.0 / 3.0 * PI;
	trq1.x = trq.x*cos(j1) - trq.y*sin(j1);
	trq1.y = trq.x*sin(j1) + trq.y*cos(j1);
	tre.x = e.x - w.x;
	tre.y = e.y - w.y;

	double anglecostrq, anglecosq, anglesintrq;
	anglesintrq = (trq1.x*(w.x - e.x) + trq1.y*(w.y - e.y)) / sqrt(trq1.x*trq1.x + trq1.y*trq1.y) / sqrt(pow(w.x - e.x, 2) + pow(w.y - e.y, 2)) - 0.5;
	anglecostrq = (trq1.x*tre.x + trq1.y*tre.y) / sqrt(trq1.x*trq1.x + trq1.y*trq1.y) / sqrt(tre.x*tre.x + tre.y * tre.y);
	anglecosq = ((q.x - w.x)*tre.x + (q.y - w.y)*tre.y) / sqrt((q.x - w.x)*(q.x - w.x) + (q.y - w.y)*(q.y - w.y)) / sqrt(tre.x*tre.x + tre.y * tre.y);

	if (anglesintrq>0 && anglecostrq < anglecosq)
	{
		tre1.x = tre.x*cos(j1) + tre.y*sin(j1);
		tre1.y = -tre.x*sin(j1) + tre.y*cos(j1);
	}
	else
	{
		trq1.x = trq.x*cos(j1) + trq.y*sin(j1);
		trq1.y = -trq.x*sin(j1) + trq.y*cos(j1);
		tre1.x = tre.x*cos(j1) - tre.y*sin(j1);
		tre1.y = tre.x*sin(j1) + tre.y*cos(j1);

	}
	trq.x = trq1.x + w.x; trq.y = trq1.y + w.y;
	tre.x = tre1.x + w.x; tre.y = tre1.y + w.y;
	POINT st=getcross(trq,e,tre,q);
	double opti;

	st.index = i;
	stainer.x = st.x; 
	stainer.y = st.y;
	cout<<" steniar point "<<st.x<<" "<<st.y<<endl;

	opti = sqrt(pow(st.x - q.x, 2) + pow(st.y - q.y, 2)) + sqrt(pow(st.x - w.x, 2) + pow(st.y - w.y, 2)) + sqrt(pow(st.x - e.x, 2) + pow(st.y - e.y, 2));
	cout<<"smt opti :"<<opti<<endl;

	double tempopti = sqrt(pow(w.x - q.x, 2) + pow(w.y - q.y, 2)) + sqrt(pow(w.x - e.x, 2) + pow(w.y - e.y, 2));
	cout<<"mst opti :"<<tempopti<<endl;

	opti = opti / tempopti;
	return opti; 
}
//stainer trique exist.ab  index-corresponded
//exist index->triangles->exist
bool  doubleadj(int a, int b, int c)	//判断是否有两条边在mst中 
{
	if (a == -2 && b == -2)
	{
		return 0;
	}
	else if (a == -2 && c == -2)
	{	
		return 0;
	}
	else if (b == -2 && c == -2)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}
double modifystainer(vector<ACS>& exist)	//优化图G 
{	
	POINT stainer;

	typedef struct cmpsmt
	{
		bool operator()(triasmt a, triasmt b)
		{
			return a.smt < b.smt;
		}
	};

	priority_queue <triasmt, vector<triasmt>, cmpsmt> trique;	//存储待处理三角形

	int j;
	for (int i = 0; i < (int)exist.size(); i++)
	{
		j = exist[i].index;
		if (!doubleadj(exist[i].ab_index, exist[i].bc_index, exist[i].ac_index))
		{
			continue;
		}
		triasmt t; 
		t.index = i; 
		t.smt = stainerpoint(triangles[j], stainer, i);
		t.stainer = stainer;
		if (t.smt < 1) 
		{
			cout<<"---->smt/mst："<<t.smt<<endl;
			trique.push(t);
		}
	}
	cout<<endl;
	j = 0;
	EDGE adj1, adj2;
	smt_cost = mst_cost;
	while (trique.size() != 0)
	{
		double lzy_cut=0;
		double cut = 0;
		triasmt t = trique.top(); 
		trique.pop();
		j = t.index;
		if (exist[j].index == -1)			//被破坏三角形 
			continue;
		if (exist[j].ab_index != -2)
		{
			if (exist[j].bc_index != -2)
			{
				adj1.index = exist[j].ab_index;

				adj1.a.index = triangles[exist[j].index].a.index;
				adj1.a.x = triangles[exist[j].index].a.x;
				adj1.a.y = triangles[exist[j].index].a.y;

				adj1.b.index = triangles[exist[j].index].b.index;
				adj1.b.x = triangles[exist[j].index].b.x;
				adj1.b.y = triangles[exist[j].index].b.y;

				adj2.index = exist[j].bc_index;

				adj2.a.index = triangles[exist[j].index].b.index;
				adj2.a.x = triangles[exist[j].index].b.x;
				adj2.a.y = triangles[exist[j].index].b.y;

				adj2.b.index = triangles[exist[j].index].c.index;
				adj2.b.x = triangles[exist[j].index].c.x;
				adj2.b.y = triangles[exist[j].index].c.y;
			}
			else
			{
				adj1.index = exist[j].ab_index;

				adj1.a.index = triangles[exist[j].index].a.index;
				adj1.a.x = triangles[exist[j].index].a.x;
				adj1.a.y = triangles[exist[j].index].a.y;

				adj1.b.index = triangles[exist[j].index].b.index;
				adj1.b.x = triangles[exist[j].index].b.x;
				adj1.b.y = triangles[exist[j].index].b.y;

				adj2.index = exist[j].ac_index;

				adj2.a.index = triangles[exist[j].index].a.index;
				adj2.a.x = triangles[exist[j].index].a.x;
				adj2.a.y = triangles[exist[j].index].a.y;

				adj2.b.index = triangles[exist[j].index].c.index;
				adj2.b.x = triangles[exist[j].index].c.x;
				adj2.b.y = triangles[exist[j].index].c.y;
			}
		}
		else
		{
			adj1.index = exist[j].bc_index;

			adj1.a.index = triangles[exist[j].index].c.index;
			adj1.a.x = triangles[exist[j].index].c.x;
			adj1.a.y = triangles[exist[j].index].c.y;

			adj1.b.index = triangles[exist[j].index].b.index;
			adj1.b.x = triangles[exist[j].index].b.x;
			adj1.b.y = triangles[exist[j].index].b.y;

			adj2.index = exist[j].ac_index;

			adj2.a.index = triangles[exist[j].index].a.index;
			adj2.a.x = triangles[exist[j].index].a.x;
			adj2.a.y = triangles[exist[j].index].a.y;

			adj2.b.index = triangles[exist[j].index].c.index;
			adj2.b.x = triangles[exist[j].index].c.x;
			adj2.b.y = triangles[exist[j].index].c.y;
		}
		TRIANGLE k = triangles[exist[j].index];

		double cut_mst=0;		//cut_smt<0

		if (G[k.a.index][k.b.index] > 0) 
		{
			cut -= G[k.a.index][k.b.index];
			cut_mst-=G[k.a.index][k.b.index];
		}
		G[k.a.index][k.b.index] = -1; 
		G[k.b.index][k.a.index] = -1;

		if (G[k.b.index][k.c.index] > 0) 
		{
			cut -= G[k.b.index][k.c.index];
			cut_mst-=G[k.b.index][k.c.index];
		}
		G[k.b.index][k.c.index] = -1;
		G[k.c.index][k.b.index] = -1;

		if (G[k.c.index][k.a.index] > 0)
		{
			cut -= G[k.c.index][k.a.index];
			cut_mst-=G[k.c.index][k.a.index];
		}
		G[k.c.index][k.a.index] = -1; 
		G[k.a.index][k.c.index] = -1;

//		cout<<"所减去的MST边的值："<<cut_mst<<endl;

		double increase_smt=0;

		//size
		for (int i = 0; i <= n; i++)
			G[i][n] = -1;
		for (int i = 0; i <= n; i++)
			G[n][i] = -1;

		if (G[adj1.a.index][n] < 0)
		{
			G[adj1.a.index][n] = sqrt(pow(adj1.a.x - t.stainer.x, 2) + pow(adj1.a.y - t.stainer.y, 2));
			G[n][adj1.a.index] = G[adj1.a.index][n]; 
			cut += G[adj1.a.index][n];
			increase_smt+=G[adj1.a.index][n];
		}

		if (G[adj1.b.index][n] < 0)
		{
			G[adj1.b.index][n] = sqrt(pow(adj1.b.x - t.stainer.x, 2) + pow(adj1.b.y - t.stainer.y, 2));
			G[n][adj1.b.index] = G[adj1.b.index][n];
			cut += G[adj1.b.index][n];
			increase_smt+=G[adj1.b.index][n];
		}

		if (G[adj2.a.index][n] < 0)
		{
			G[adj2.a.index][n] = sqrt(pow(adj2.a.x - t.stainer.x, 2) + pow(adj2.a.y - t.stainer.y, 2));
			G[n][adj2.a.index] = G[adj2.a.index][n]; 
			cut += G[adj2.a.index][n];
			increase_smt+=G[adj2.a.index][n];
		}

		if (G[adj2.b.index][n] < 0)
		{
			G[adj2.b.index][n] = sqrt(pow(adj2.b.x - t.stainer.x, 2) + pow(adj2.b.y - t.stainer.y, 2));
			G[n][adj2.b.index] = G[adj2.b.index][n]; 
			cut += G[adj2.b.index][n];
			increase_smt+=G[adj2.b.index][n];
		}

//		cout<<"所增加的SMT边的值："<<increase_smt<<endl;

		points.push_back(t.stainer);			//加入点集 
		int tcv = adj1.index;
		if (adj1.index >= 0) 
		{
			exist[triangles[tcv].MST_index].index = -1;
		}

		tcv = adj2.index;
		if (adj2.index >= 0) 
		{
			exist[triangles[tcv].MST_index].index = -1;
		}

		n++;
		
		lzy_cut=increase_smt+cut_mst;
		cout<<"SMT局部优化："<<lzy_cut<<endl;
//		cout<<"sssssssssssssssssssssss  "<<cut<<endl;
		smt_cost += lzy_cut;
	}
	cout << "\nsmt:\nsmt cost:" << smt_cost << endl;
	return 0;
}

int main(void)
{
	double smt;
	Input();
	Delaunay();
	MST();
	adjcent();
	vector<ACS> exisarc(acs); 
	modifystainer(exisarc);
	cout << "\nsmt/mst ："<<smt_cost / mst_cost<<"\n\n";
	system("pause");
	return 0;
}