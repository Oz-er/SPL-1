#include<bits/stdc++.h>
using namespace std;
typedef long long ll;


struct Cereal{
    string name;
    ll calories;
    ll protein;
    ll fat;


    double mCalories;
    double mProtein;
    double mFat;
};



void calculateMarginal(vector<Cereal>&cereals, map<ll,ll>&mCal,map<ll,ll>&mProtein,map<ll,ll>&mFat){

    ll n = cereals.size();

    double total  =0;

    for(ll i=0;i<n;i++){
        cereals[i].mCalories = (double)mCal[cereals[i].calories]/n ;
        cereals[i].mProtein = (double)mProtein[cereals[i].protein]/n ;
        cereals[i].mFat = (double)mFat[cereals[i].fat]/n ;
    }


    for(auto& p :mCal){
        cout<<"calories "<<p.first<<" : "<<(double)p.second/n<<endl;
    }
    for(auto& p :mProtein){
        cout<<"protein "<<p.first<<" : "<<(double)p.second/n<<endl;
    }
    for(auto& p :mFat){
        cout<<"fat "<<p.first<<" : "<<(double)p.second/n<<endl;
    }



}


void inputFromCSV(vector<Cereal>&cereals, map<ll,ll>&mCal,map<ll,ll>&mProtein,map<ll,ll>&mFat , map<ll,map<ll,ll>> &jCalcProtein,map<ll,map<ll,ll>> &jProteinFat,map<ll,map<ll,ll>> &jCalcFat){

    ifstream file1;

    string line = "";


    file1.open("cerealinfo.csv");
    if(!file1.is_open()){
    cerr << "Error: Cannot open cerealinfo.csv\n";
    exit(1);
    }


    getline(file1, line);




    while(getline(file1,line)){

        stringstream ss(line);
        string name,temp;
        Cereal c;


        getline(ss,name,',');
        c.name=name;

        getline(ss,temp,',');
        c.calories=stoll(temp);
        mCal[c.calories]++;
        
  
        getline(ss,temp,',');
        c.protein=stoll(temp);
        mProtein[c.protein]++;

        getline(ss,temp,',');
        c.fat=stoll(temp);
        mFat[c.fat]++;


        jCalcProtein[c.calories][c.protein]++;
        jProteinFat[c.protein][c.fat]++;
        jCalcFat[c.calories][c.fat]++;

        cereals.push_back(c);

    }

    file1.close();
}


void calculateConditional(vector<Cereal>&cereals, map<ll,ll>&mCal,map<ll,ll>&mProtein,map<ll,ll>&mFat,map<ll,map<ll,ll>> &jCalProtein,map<ll,map<ll,ll>> &jProteinFat,map<ll,map<ll,ll>> &jCalcFat , ll x, ll y, ll mode){

    ll num;
    ll denum;
    
    if(mode ==1){
        num = jCalProtein[x][y];
        denum=mProtein[y];
    }
    else if(mode ==2){
        num = jProteinFat[x][y];
        denum=mFat[y];
    }
    else if(mode ==3){
        num = jCalcFat[x][y];
        denum=mFat[y];
    }
    else if(mode ==4){
        num = jCalProtein[x][y];
        denum=mProtein[x];
    }
    else if(mode ==5){
        num = jProteinFat[x][y];
        denum=mFat[x];
    }
    else if(mode ==6){
        num = jCalcFat[x][y];
        denum=mFat[x];
    }

    else {
        cout << "Invalid mode" << endl;
        return;
    }

    if(denum==0){
        cout<<"NOT ENOUGH DATA"<<endl;
    }

    double prob = (double)num / (double)denum;
    cout <<"Conditional probability P(A|B) ="<<prob<<endl;
    
}



int main(){

    map<ll,ll> mCal,mProtein,mFat;
    vector<Cereal> cereals;
    map<ll,map<ll,ll>> jCalProtein,jProteinFat,jCalFat;





    inputFromCSV(cereals,mCal,mProtein,mFat,jCalProtein,jProteinFat,jCalFat);



    cout<<"1.Conditional Probability"<<endl;
    cout<<"2.Marginal Probability"<<endl;

    ll choice;
    cin>>choice;

    if(choice==1){
        calculateMarginal(cereals,mCal,mProtein,mFat);
    }

    else if(choice ==2){

    cout<<"INPUT MODE : "<<endl;
    cout<<"1.Calorie|Protein"<<endl;
    cout<<"2.Protein|Fat"<<endl;
    cout<<"3.Calories|Fat"<<endl;
    cout<<"4.Protein|Calories"<<endl;
    cout<<"5.Fat|Protein"<<endl;
    cout<<"6.Fat|Calories"<<endl;

    ll c;
    cin>>c;

    cout<<"probability of A|B"<<endl;
    cout<<"A : ";
    ll a;
    cin>>a;
    cout<<"B : ";
    ll b;
    cin>>b;


    calculateConditional(cereals,mCal,mProtein,mFat,jCalProtein,jProteinFat,jCalFat,a,b,c);
    }


    else{
        cout<<"Invalid input"<<endl;
    }



}