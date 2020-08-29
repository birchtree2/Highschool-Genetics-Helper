#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cctype>
#include<map>
#include<cmath>
#include<vector>
#include "fraction.h"
int genome_cnt(std::string s) { //染色体组数
	for(int i=1; i<(int)s.length(); i++) if(tolower(s[i])!=tolower(s[0])) return i;
	return (int)s.length();
}
int gene_cnt(std::string s) { //等位基因对数
	return (int)s.length()/genome_cnt(s);
}
void _gametes(std::vector<std::string>& ans,std::string fa,std::string s,int deep) {
	if(deep==(int)fa.length()) {
		ans.push_back(s);
		return;
	}
	int gn=genome_cnt(fa);
	for(int i=deep; i<=deep+gn-1; i++) _gametes(ans,fa,s+fa[i],deep+gn);
}
std::vector<std::string> gametes(std::string s) { //生成配子
	std::vector<std::string>ans;
	_gametes(ans,s,"",0);
	return ans;
}
std::string forma(std::string s) { //按照基因的写法排序
	std::sort(s.begin(),s.end(),
	[](char p,char q)->bool {
		if(tolower(p)==tolower(q)) return p<q;
		else return tolower(p)<tolower(q);
	});
	return s;
}
std::string phenotype(std::string s) {
	std::string ans="";
	int gn=genome_cnt(s);
	for(int i=0; i<(int)s.length(); i+=gn) {
		bool is_recess=1;//是否表现为隐性
		for(int j=0; j<gn; j++) if(s[i+j]!=tolower(s[i])) is_recess=0;
		if(is_recess) ans+=tolower(s[i]);
		else ans+=toupper(s[i]);
	}
	return ans;
}
struct Genepool {
	std::map<std::string, Fraction> pr;//每种基因型的概率
	void print_ratio(std::map<std::string,Fraction> ma) {
		for(auto it=ma.begin(); it!=ma.end(); it++) {
			std::cout<<it->first;
			auto p=it;
			p++;
			if(p==ma.end()) std::cout<<"=";
			else std::cout<<":";
		}
		int lcmdown=ma.begin()->second.getDenominator();//通分分母
		for(auto it=(++ma.begin()); it!=ma.end(); it++) {
			int down = ma.begin()->second.getDenominator();
			lcmdown = lcm(lcmdown, down);
		}
		for(auto it=ma.begin(); it!=ma.end(); it++) {
			std::cout<<lcmdown/it->second.getDenominator()*it->second.getNumerator();
			auto p=it;
			p++;
			if(p!=ma.end())std::cout<<":";
		}
	}
	void analyze() {
		std::map<std::string,Fraction >pheno_pr;
		for(auto s : pr) {
			pheno_pr[phenotype(s.first)] += s.second;
		}
		std::cout<<"表现型比例\n ";
		print_ratio(pheno_pr);
		std::cout << "\n";
		std::cout << "基因型比例:\n";
		print_ratio(pr);
		std::cout << "\n";

	}
};
Genepool cross(std::string fa,std::string mo) {
	std::vector<std::string>gfa,gmo;
	Genepool ans;
	gfa=gametes(fa);
	gmo=gametes(mo);
	int sz = (int)gfa.size() * (int)gmo.size();
	for(std::string p : gfa) {
		for(std::string q : gmo) {
			ans.pr[forma(p + q)] += Fraction(1, sz);
//        	std::cout<<Fraction(1,sz)<<endl;
//        	std::cout<<ans.pr[forma(p + q)]<<endl;
		}
	}
	return ans;
}


Genepool breed(Genepool orig,int type) {
	Genepool ans;
	Genepool tmp;
	if(type==1) { //自交
		for(auto s : orig.pr) {
			tmp=cross(s.first,s.first);
			for(auto t : tmp.pr) {
				ans.pr[t.first] += t.second * s.second;
			}
		}
	} else if(type==2) { //自由交配
		for(auto p : orig.pr) {
			for(auto q : orig.pr) {
				tmp=cross(p.first,q.first);
				for(auto t : tmp.pr) {
					ans.pr[t.first] += t.second * p.second * q.second;
				}
			}
		}
	}
	return ans;
}
int main() {
//	Fraction a=Fraction(0,0)+Fraction(1,4);
//	a-=b;
//	std::cout<<a;
	std::string s1,s2;
	std::vector<std::string>in,out;
	while(std::cin>>s1>>s2) {
		cross(s1,s2).analyze();
//		in.push_back(s1);
//		in.push_back(s2);
//		out=breed(in,1);
//		analyze(out);
	}
	system("pause");
}

