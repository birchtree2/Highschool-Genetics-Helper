#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cctype>
#include<map> 
#include<cmath> 
int genome_cnt(std::string s){//染色体组数 
	for(int i=1;i<(int)s.length();i++) if(tolower(s[i])!=tolower(s[0])) return i;
	return (int)s.length();
}
int gene_cnt(std::string s){//等位基因对数 
	return (int)s.length()/genome_cnt(s);
}
void _gametes(std::vector<std::string>& ans,std::string fa,std::string s,int deep){
	if(deep==(int)fa.length()){
		ans.push_back(s);
		return; 
	}
	int gn=genome_cnt(fa);
	for(int i=deep;i<=deep+gn-1;i++) _gametes(ans,fa,s+fa[i],deep+gn);
}
std::vector<std::string> gametes(std::string s){//生成配子 
	std::vector<std::string>ans;
	_gametes(ans,s,"",0);
	return ans;
}
std::string forma(std::string s){//按照基因的写法排序 
	std::sort(s.begin(),s.end(),
	[](char p,char q)->bool{
		if(tolower(p)==tolower(q)) return p<q;
		else return tolower(p)<tolower(q);
	});
	return s;
} 
std::vector<std::string> cross(std::string fa,std::string mo){
	std::vector<std::string>gfa,gmo,ans;
	gfa=gametes(fa);
	gmo=gametes(mo);
	for(std::string p : gfa){
		for(std::string q : gmo) ans.push_back(forma(p+q));
	}
	return ans;
}
std::string phenotype(std::string s){
	std::string ans="";
	int gn=genome_cnt(s);
	for(int i=0;i<(int)s.length();i+=gn){
		bool is_recess=1;//是否表现为隐性 
		for(int j=0;j<gn;j++) if(s[i+j]!=tolower(s[i])) is_recess=0;
		if(is_recess) ans+=tolower(s[i]);
		else ans+=toupper(s[i]);
	}
	return ans;
} 
void print_ratio(std::map<std::string,int> ma){
	for(auto it=ma.begin();it!=ma.end();it++){
		std::cout<<it->first;	
		auto p=it;
		p++;
		if(p==ma.end()) std::cout<<"=";
		else std::cout<<":"; 
	} 
 	int g=ma.begin()->second;
 	for(auto it=(++ma.begin());it!=ma.end();it++) g=std::__gcd(g,it->second);
	for(auto it=ma.begin();it!=ma.end();it++){
		std::cout<<it->second/g;	
		auto p=it;
		p++;
		if(p!=ma.end())std::cout<<":"; 
	} 
} 
void analyze(std::vector<std::string> vec){
	std::map<std::string,int>gene_num;
	std::map<std::string,int>pheno_num;
	for(std::string s : vec){
		gene_num[s]++; 
		pheno_num[phenotype(s)]++;
	} 
	std::cout<<"表现型比例:\n";
	print_ratio(pheno_num);
	std::cout<<"\n";
	std::cout<<"基因型比例:\n";
	print_ratio(gene_num);
	std::cout<<"\n"; 
	 
}
std::vector<std::string> breed(std::vector<std::string> vec,int type){
	vector<std::string>ans;
	vector<std::string>tmp;
	if(type==1){//自交 
		for(auto s : vec){
			tmp=cross(s,s);
			for(auto t : tmp) ans.push_back(t);
		}  
	}else if(type==2){//自由交配 
		for(auto p : vec){
			 for(auto q : vec){
			 	tmp=cross(p,q);
			 	for(auto t : tmp) ans.push_back(t); 
			 }
		} 
	}
	return ans;
} 
int main(){
	std::string s1,s2;
	while(std::cin>>s1>>s2){
		analyze(cross(s1,s2));
	} 
}

