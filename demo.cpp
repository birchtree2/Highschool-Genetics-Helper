#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cctype>
#include<map>
#include<vector>
#include "fraction.h"
int genome_cnt(std::string s) { //Ⱦɫ������
	for(int i=1; i<(int)s.length(); i++) if(tolower(s[i])!=tolower(s[0])) return i;
	return (int)s.length();
}
int gene_cnt(std::string s) { //��λ�������
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
std::vector<std::string> gametes(std::string s) { //��������
	std::vector<std::string>ans;
	_gametes(ans,s,"",0);
	return ans;
}
std::string forma(std::string s) { //���ջ����д������
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
		bool is_recess=1;//�Ƿ����Ϊ����
		for(int j=0; j<gn; j++) if(s[i+j]!=tolower(s[i])) is_recess=0;
		if(is_recess) ans+=tolower(s[i]);
		else ans+=toupper(s[i]);
	}
	return ans;
}
struct Genepool {
	std::map<std::string, Fraction> pr;//ÿ�ֻ����͵ĸ���
	Genepool(){}
	Genepool(std::string s){
		pr[s]=1;
	}
	Genepool(std::vector<std::string> vec){
		for(auto s : vec) pr[s]+=Fraction(1,vec.size());
	}
	void print_ratio(std::map<std::string,Fraction> ma) {
		for(auto it=ma.begin(); it!=ma.end(); it++) {
			std::cout<<it->first;
			auto p=it;
			p++;
			if(p==ma.end()) std::cout<<"=";
			else std::cout<<":";
		}
		int lcmdown=ma.begin()->second.getDenominator();//ͨ�ַ�ĸ
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
		std::cout<<"�����ͱ���\n ";
		print_ratio(pheno_pr);
		std::cout << "\n";
		std::cout << "�����ͱ���:\n";
		print_ratio(pr);
		std::cout << "\n";

	}
	void clear(){
		pr.clear();
	}
};
Genepool cross(std::string fa,std::string mo) {
	std::vector<std::string>gfa,gmo,ans;
	gfa=gametes(fa);
	gmo=gametes(mo);
	// int sz = (int)gfa.size() * (int)gmo.size();
	for(std::string p : gfa) {
		for(std::string q : gmo) {
			ans.push_back(forma(p+q));
		}
	}
	return Genepool(ans);
}


Genepool breed(Genepool parent,int generation,int type) {
	Genepool ans;
	Genepool tmp;
	for(int tim=1;tim<=generation;tim++){
		ans.clear();
		if(type==1) { //�Խ�
			for(auto s : parent.pr) {
				tmp=cross(s.first,s.first);
				for(auto t : tmp.pr) {
					ans.pr[t.first] += t.second * s.second;
				}
			}
		} else if(type==2) { //���ɽ���
			for(auto p : parent.pr) {
				for(auto q : parent.pr) {
					tmp=cross(p.first,q.first);
					for(auto t : tmp.pr) {
						ans.pr[t.first] += t.second * p.second * q.second;
					}
				}
			}
		}
		parent=ans;
	}
	return ans;
}
int main() {
//	Fraction a=Fraction(0,0)+Fraction(1,4);
//	a-=b;
//	std::cout<<a;
	std::string s1,s2;
	std::vector<std::string>in,out;
	Genepool ans;
	while(std::cin>>s1) {
// cross(s1,s2).analyze();
		in.push_back(s1);
		ans=Genepool(in);
		ans=breed(ans,2,1);
		ans.analyze();
//		out=breed(in,1);
//		analyze(out);
	}
	system("pause");
}

