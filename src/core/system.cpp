/*************************************************************************************/
/*      Copyright 2009 Barcelona Supercomputing Center                               */
/*                                                                                   */
/*      This file is part of the NANOS++ library.                                    */
/*                                                                                   */
/*      NANOS++ is free software: you can redistribute it and/or modify              */
/*      it under the terms of the GNU Lesser General Public License as published by  */
/*      the Free Software Foundation, either version 3 of the License, or            */
/*      (at your option) any later version.                                          */
/*                                                                                   */
/*      NANOS++ is distributed in the hope that it will be useful,                   */
/*      but WITHOUT ANY WARRANTY; without even the implied warranty of               */
/*      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                */
/*      GNU Lesser General Public License for more details.                          */
/*                                                                                   */
/*      You should have received a copy of the GNU Lesser General Public License     */
/*      along with NANOS++.  If not, see <http://www.gnu.org/licenses/>.             */
/*************************************************************************************/

#include "system.hpp"
#include "config.hpp"
#include "plugin.hpp"
#include "schedule.hpp"
#include "barrier.hpp"
#include "nanos-int.h"
#include "copydata.hpp"
#include "os.hpp"
#include "basethread.hpp"
#include "malign.hpp"
#include "processingelement.hpp"
#include "allocator.hpp"
#include "debug.hpp"
#include "dlb.hpp"
#include <string.h>
#include <set>
#include <climits>

#ifdef SPU_DEV
#include "spuprocessor.hpp"
#endif

#ifdef GPU_DEV
#include "gpuprocessor_decl.hpp"
#endif

#ifdef OpenCL_DEV
#include "openclprocessor.hpp"
#endif

#ifdef PAPI
#include <iostream>
#include <fstream>
#include <limits>
#include <iomanip>

#endif

using namespace nanos;

System nanos::sys;

#ifdef PAPI
void System::addWdDesToSet(std::string wddescription){

//<>std::cout << "adding to set "<< wddescription <<std::endl;
	std::string aux=wddescription;
	
	labelconf newl;// = new labelconf();
	
	newl.slabel=wddescription;
	newl.mode=0;
	newl.numconsec=0;	
	{
	//->LockBlock lock(_hashlock);//new
	_wdDescriptionMmapConf.insert(std::pair<std::string,labelconf>(wddescription,newl));
	//LockBlock lock2(_hashlock);//newer
  
	_wdDescriptionSet.insert(wddescription);
	}
}
//////////////////////////////////Mmap wdDescriptionsetconf/////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
void System::insertTaskToMmapconf(std::string wddescription){
{	
//->LockBlock lock(_hashlock2);//new
labelconf lc;
lc.slabel=wddescription; lc.mode=0; lc.numconsec=0; lc.prefetcherconf=0; lc.numchecks=0;
////std::cout<<"--->Inserted to wdDescriptionMmap.conf: " << lc.prefetcherconf <<" "<< lc.numconsec<<  std::endl;
_wdDescriptionMmapConf.insert ( std::pair<std::string,labelconf>(wddescription,lc) );
	clearBLock(wddescription);
	_wdDescriptionSet.insert(wddescription);//new
}
}

//checks whether a task of a given label already exists in the set, (already has had an instance)
//true exists, false not
bool System::checkTaskToMmapConf(std::string wddescription){
	//std::set<std::string>::iterator it= _wdDescriptionSet.find ( wddescription  );
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	if(it == _wdDescriptionMmapConf.end()){
		return false;
	}else{
		return true;
	}
}
//increase num of consecutive times that task type has been seen
void System::increaseMmapConfCont(std::string wddescription){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	//if(it->second.mode==0)
		{
		//->LockBlock lock(_hashlock3);   ///whatch out!!!!
		//<>std::cout << "I'm thread: " << myThread->getId() << " numconsec was: " << it->second.numconsec <<std::endl;
		it->second.numconsec=(it->second.numconsec+1);//%2;
		}
	//else if (it->second.mode==1)
	//	it->second.numconsec=(it->second.numconsec+1)%10;
}
//----------------bLock----------------
void System::tangleBLock( std::string wddescription, int thread ){
std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
for (std::multimap<std::string,labelconf>::iterator it2=_wdDescriptionMmapConf.begin(); it2!=_wdDescriptionMmapConf.end(); ++it2){
	if(it2->second.slabel.compare(it->second.slabel)!=0){untangleBLock(wddescription, thread);}}
it->second.bLock[thread]=true;
}
void System::untangleBLock( std::string wddescription, int thread ){
std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
it->second.bLock[thread]=false;
}
void System::clearBLock(  std::string wddescription ){//should be called at the begining
std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
for(int i=0;i<sys.getNumThreads();i++)
	it->second.bLock[i]=false;
}
bool System::checkPosBLock (std::string wddescription, int thread){
std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
return it->second.bLock[thread];
}
bool System::checkBLock (std::string wddescription){//return true if all positions are true
std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription); bool all=true;
for(int i=0;i<sys.getNumThreads() && all==true ;i++){
	if(it->second.bLock[i]==false) {all=false;}} return all;
}
//-----------pref type task consec
//resets the numconsec counter. starting at one!!
void System::resetMmapConfCont(std::string wddescription){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
		it->second.numconsec=0;
		////std::cout<<"configuration reseted confcont " << std::endl;
}
//make sure you use this function when the wddescription exists
int System::checkNumConsecMmapConf(std::string wddescription){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	////std::cout<<"numconsec:  " << it->second.numconsec << std::endl;
	return it->second.numconsec;
}
//------------pref conf mode
//make sure you use this function when the wddescription exists
void System::changeModeMmapConf(std::string wddescription, int mode){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	it->second.mode=mode;
}
int System::checkModeMmapConf(std::string wddescription){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	return it->second.mode;
}
//-------pref conf number
void System::setMmapPrefConf(int preconf, std::string wddescription){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	it->second.prefetcherconf=preconf;
}
int System::checkMmapPrefConf(std::string wddescription){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	return it->second.prefetcherconf;
}
void System::increaseMmapPrefConf(std::string wddescription){//to start with we'll have 5 confs from 0 to 4
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	if(it->second.prefetcherconf==5)
		it->second.prefetcherconf=0;
	else
		it->second.prefetcherconf=(it->second.prefetcherconf+1);
	
	////std::cout<<"numconfpref:  " << it->second.prefetcherconf << std::endl;
}
void System::resetMmapPrefConf(std::string wddescription){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	it->second.prefetcherconf=0;
	////std::cout<<"conf reseted:  " << it->second.prefetcherconf << std::endl;
}

std::string System::translateNumberToReal(int reg){

std::string regc;
if(reg==0)//disabled 001
  regc="1";
else if(reg==1)//default
  regc="0";
else if(reg==2)//shallowest
  regc="2";
else if(reg==3)//deepest
  regc="7";
else if(reg==4)//prefetch on stores
  regc="8";
else if(reg==5)//most aggressive
  regc="1f";

//std::cout<< "translating...: " << regc << std::endl;	
return regc;
}
/////
//-------num checkings
int System::checkMmapCountCent(std::string wddescription ){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	return it->second.numchecks;
}
void System::increaseMmapCountcent(std::string wddescription){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	it->second.numchecks++;
}
void System::resetMmapCountCent(std::string wddescription){
	std::multimap<std::string,labelconf>::iterator it = _wdDescriptionMmapConf.find(wddescription);
	it->second.numchecks=0;
}

//-------------
int System::choseBestPrefConf(std::string wddescription){
std::set<std::string>::iterator it;
it=_wdDescriptionSet.find( wddescription);
std::pair <std::multimap<std::string ,nanos::WorkDescriptor::wdinfo >::reverse_iterator, std::multimap<std::string,nanos::WorkDescriptor::wdinfo>::reverse_iterator> ret;
    	ret = _hashSysStatistics.equal_range(*it);

	//<>std::cout<< "    ->length hash: " <<_hashSysStatistics.count(wddescription) << std::endl;

			//labelstats aux;//=new labelstats();
			//aux.slabel= *it;
	int i=0; //gather last 10 ipc from a given type of task
	float bestipc=0.0;
	//float medcurripc=0.0;
	int  bestpreconf=0;
/*	for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::reverse_iterator it2=ret.second; it2!=ret.first and i<80; ++it2){	
		std::cout<<"conf ZZ: "<< it2->second.prefetcherconf <<" "<< it2->second.ipc <<std::endl;
		i++;
	}*/
i=0;
float ipcs[6];
int ipcsfound[6];

		for(int k=0;k<6;k++){
			ipcs[k]=0.0;
			ipcsfound[k]=0;
		}

	
	//<>std::cout <<"before  done the loop " << std::endl;
	bool getout=false;
// i<6*sys.getNumPSearch()//condition beforetherewas getout var
	for (std::multimap<std::string,nanos::WorkDescriptor::wdinfo>::reverse_iterator it2=ret.second; !getout && it2!=ret.first ; ++it2){		
		//<>std::cout << "loop: conf:" << it2->second.prefetcherconf <<" i "<<i<< std::endl;

		if(ipcsfound[it2->second.prefetcherconf]<sys.getNumPSearch()){
			ipcs[it2->second.prefetcherconf]+=it2->second.ipc;
			//std::cout<<"Conf: "<< it2->second.prefetcherconf << " IPC: " << it2->second.ipc << std::endl;
			ipcsfound[it2->second.prefetcherconf]+=1;	
		}	
		i++;
		bool alarm=false;
		for(int l=0; l<6 && !alarm ;l++){
			if(ipcsfound[l]<sys.getNumPSearch()){
				alarm=true;
			}
		}	
		if ( !alarm ){
				//<>std::cout<< "!!!!!!!!ALARM in false" << std::endl;
				//<>for(int l=0; l<6; l++){
				//<>	std::cout << "value alarm i: " << ipcsfound[l] << " "<< l << std::endl;}
			 getout=true;}
	}
//std::cout <<" done the loop " << std::endl;
		//bestipc=ipcs[0]/(float)sys.getNumPSearch();
		if (ipcsfound[0]>0){
			bestipc=ipcs[0]/(float)ipcsfound[0];}
		
//std::cout <<" done the loop 2" << std::endl;
		bestpreconf=0;
		//for(int k=1;k<6;k++){
		int k=1;
		while( k<6){
			if(k==1){k=2;}
			else if(k==2){k=1;}
	
			if(ipcsfound[k]>0){
				//ipcs[k]/=(float)sys.getNumPSearch();
				ipcs[k]/=(float)ipcsfound[k];
			}
			if( (ipcs[k] - (_kbPowerBW * 0.01)) > bestipc){
				bestipc=ipcs[k];
				bestpreconf=k;}


		if(k==2){k=1;}
		else if(k==1){k=2;}
		k++;
		}
//<>std::cout <<"                   choosen conf " << bestpreconf <<std::endl;
		/*if(bestpreconf==0)
			bestpreconf=1;
		else if (bestpreconf==1)
			bestpreconf=0;
		else if (bestpreconf==2)
			bestpreconf=7;
		else if (bestpreconf==3)
			bestpreconf=8;*/
		//else if (bestpreconf==4)
		//	bestpreconf=16;
		//std::cout<<"best conf: "<< bestpreconf <<std::endl;
		return bestpreconf;
}
/*int System::choseBestPrefConf(std::string wddescription){
std::set<std::string>::iterator it;
it=_wdDescriptionSet.find( wddescription);
std::pair <std::multimap<std::string ,nanos::WorkDescriptor::wdinfo >::reverse_iterator, std::multimap<std::string,nanos::WorkDescriptor::wdinfo>::reverse_iterator> ret;
    	ret = _hashSysStatistics.equal_range(*it);
			//labelstats aux;//=new labelstats();
			//aux.slabel= *it;
	int i=0,j=0; //gather last 10 ipc from a given type of task
	float bestipc=0.0;
	float medcurripc=0.0;
	int antpreconf=4, currpreconf=0, bestpreconf=0;
	
	for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::reverse_iterator it2=ret.second; it2!=ret.first and i<10; ++it2){	
			
		////std::cout<<"conf ZZ: "<< it2->second.prefetcherconf <<" "<< it2->second.ipc <<std::endl;
		i++;
	}

i=0;j=0;
	for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::reverse_iterator it2=ret.second; i<10; ++it2){		
		if(i==0)	
				antpreconf=it2->second.prefetcherconf;
		currpreconf=it2->second.prefetcherconf;
		////std::cout<<"conf: "<< currpreconf <<std::endl;
		if(currpreconf==antpreconf){//j<=1) {//confs from the same label are the same
					medcurripc+=it2->second.ipc;
					j++;//should take the value of the number of times a conf is tested
					//std::cout<<"-equal-"<<std::endl;
		}else{//compute the median of the previous, update best, and start to accumulate the current
				////std::cout<<"medcurripc is: " << medcurripc <<" "<< medcurripc/float(j) <<" for conf: "<< antpreconf <<std::endl;
				medcurripc=medcurripc/float(j);
				j=1;
				if (medcurripc > bestipc) {
					bestipc=medcurripc;
					medcurripc=0.0;
					bestpreconf=currpreconf;
				}
				medcurripc=it2->second.ipc;
		}
		antpreconf=currpreconf;
		i++;
		//j++;
	}
	//j++;
	////std::cout<<"medcurripc is: " << medcurripc <<" "<< medcurripc/float(j) <<" for conf: "<< antpreconf <<std::endl;
	medcurripc=medcurripc/float(j);
	if (medcurripc > bestipc) {
					bestipc=medcurripc;
					bestpreconf=antpreconf;
	}

		////std::cout<<"BEST configuration is: " <<bestpreconf<< std::endl;
		return bestpreconf;
}
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void System::addResultToGroups(WorkDescriptor::wdinfo data){
	//->LockBlock lock(_hashlock4);
	//std::cout<< "adding to groups " << data.prefetcherconf <<std::endl;
	
	_hashSysStatistics.insert(std::make_pair(data.wdlabel, data ));
	//std::cout<< "added to groups" <<std::endl;
}

void System::addResultToStatistics(WorkDescriptor::wdinfo data){
	//->LockBlock lock(_hashlock5);
	_sysStatistics.push_back(data);
}
//prints stdout each task type statistics
void System::printGroups(){

  std::cout << "printing hash :\n";
  //std::set<std::string>::iterator it;
  std::set<std::string>::iterator it;
  for (it=_wdDescriptionSet.begin(); it!=_wdDescriptionSet.end(); ++it){

    //std::pair <std::multimap<std::string ,nanos::WorkDescriptor::wdinfo >::iterator, std::multimap <std::string,nanos::WorkDescriptor::wdinfo>::iterator> ret;
   
    std::pair <std::multimap<std::string ,nanos::WorkDescriptor::wdinfo >::iterator, std::multimap <std::string,nanos::WorkDescriptor::wdinfo>::iterator> ret;

    ret = _hashSysStatistics.equal_range(*it);
    std::cout << *it << " =>\n";
	

    for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::iterator it2=ret.first; it2!=ret.second; ++it2){
	
	 std::cout << "value of numInstance: " << _numIns  << std::endl;
	//std::cout.precision(15);
	//std::cout.precision( std::numeric_limits< floati >::digits10);
      	std::cout << " ipc " <<  std::fixed << std::setprecision(22) << it2->second.ipc;
	std::cout << " ins " << it2->second.ins;
	std::cout << " l1tcm " << it2->second.l1tcm;
	std::cout << " dins " << it2->second.dins;
	std::cout << "\n";
    }

    std::cout << '\n';
  }

}

void System::writeDump(){
	std::cout << "Dumping"<< std::endl;
	std::ofstream fileipc,fileinst;
	
	char threads[15]; sprintf(threads, "%d", _numThreads);
	char numps[15]; sprintf(numps, "%d", _numPS);
	char numpb[15]; sprintf(numpb, "%d", _numPB);
  char dot[1]; dot[0]='.'; char stripc[61]; char strins[61];
	//strcpy (str,exec);
	//strcat (str,dot);
	strcpy (stripc,threads); strcat (stripc,dot); strcat (stripc,numps); strcat (stripc,dot); strcat (stripc,numpb);
	strcat (stripc,dot); //strcat (str,"statistics");
	strcpy (strins,threads); strcat (strins,dot); strcat (strins,numps); strcat (strins,dot); strcat (strins,numpb);
	strcat (strins,dot); //strcat (str,"statistics");

  if(_kbPowerBW>0) {
                char kpbw [15];//knobpowerbandwidth
                sprintf(kpbw, "%d", _kbPowerBW);
								strcat(stripc,kpbw);strcat(stripc,dot);
								strcat(strins,kpbw);strcat(strins,dot);
  }
	if(_numIns>0) {
                char numinstance [15];
                sprintf(numinstance, "%d", _numIns);
                //std::string numinstancest=numinstance;
                //pref=pref+numinstancest;
                //pref=pref+dot;
								strcat(stripc,numinstance);strcat(stripc,dot);
								strcat(strins,numinstance);strcat(strins,dot);
  }

	strcat(stripc,"IPC.txt");
	strcat(strins,"PAPI_TOT_INS.txt");

  	fileipc.open  ((char*)stripc, std::ios::trunc );//app
	fileinst.open ((char*)strins, std::ios::trunc );

	std::set<std::string>::iterator it;
        //for each label
        for (it=_wdDescriptionSet.begin(); it!=_wdDescriptionSet.end(); ++it){

                std::pair <std::multimap<std::string ,nanos::WorkDescriptor::wdinfo >::iterator, std::multimap<std::string,nanos::WorkDescriptor::wdinfo>::iterator> ret;

		ret = _hashSysStatistics.equal_range(*it);

		//std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::iterator it3=ret.first;
		//if (it3!=NULL) myfile << getNumThreads()<<":"<< it3->second.slabel;
		
		int size= _hashSysStatistics.size();

		if(size!=0){

			bool cent=true;
			for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::iterator it2=ret.first; it2!=ret.second; ++it2){
				if(cent){ 
					fileipc << getNumThreads() << ":" <<  it2->second.wdlabel; cent=false;
					fileinst << getNumThreads() << ":" <<  it2->second.wdlabel; cent=false;
					}
	
				fileipc << ":" << it2->second.ipc;
				fileinst << ":" << it2->second.ins;

			}
		fileipc  << std::endl;
		fileinst << std::endl;
		}
	}
	//fileipc << getNumThreads()<<":ipc:"<<meantotalipc<< std::endl;
		
}

void System::writeDumpConfs(){
	std::cout << "Dumping"<< std::endl;
	std::ofstream fileipc,fileinst,fileipcstable;
	
	char threads[15]; sprintf(threads, "%d", _numThreads);
	char numps[15]; sprintf(numps, "%d", _numPS);
	char numpb[15]; sprintf(numpb, "%d", _numPB);

  char dot[1]; dot[0]='.'; char stripc[61]; char strins[61]; char stripcstable[61];
	//strcpy (str,exec);
	//strcat (str,dot);
	strcpy (stripc,threads); strcat (stripc,dot); strcat (stripc,numps); strcat (stripc,dot); strcat (stripc,numpb);
	strcat (stripc,dot); //strcat (str,"statistics");

	strcpy (strins,threads); strcat (strins,dot); strcat (strins,numps); strcat (strins,dot); strcat (strins,numpb);
	strcat (strins,dot); //strcat (str,"statistics");

	strcpy (stripcstable,threads); strcat (stripcstable,dot); strcat (stripcstable,numps); strcat (stripcstable,dot); strcat (stripcstable,numpb);
	strcat (stripcstable,dot); //strcat (str,"statistics");

	if(_kbPowerBW>0) {
                char kpbw [15];//knobpowerbandwidth
                sprintf(kpbw, "%d", _kbPowerBW);
								strcat(stripc,kpbw);strcat(stripc,dot);
								strcat(strins,kpbw);strcat(strins,dot);
								strcat(stripcstable,kpbw);strcat(stripcstable,dot);
  }

	if(_numIns>0) {
                char numinstance [15];
                sprintf(numinstance, "%d", _numIns);
                //std::string numinstancest=numinstance;
                //pref=pref+numinstancest;
                //pref=pref+dot;
								strcat(stripc,numinstance);strcat(stripc,dot);
								strcat(strins,numinstance);strcat(strins,dot);
								strcat(stripcstable,numinstance);strcat(stripcstable,dot);
  }

	strcat(stripc,"IPCconfs.txt");
	strcat(stripcstable,"IPCstable.txt");
	strcat(strins,"PAPI_TOT_INSconfs.txt");

  fileipc.open  ((char*)stripc, std::ios::trunc );//app
	fileinst.open ((char*)strins, std::ios::trunc );
	fileipcstable.open ((char*)stripcstable, std::ios::trunc );

		double meanipcconf0=0.0;
		double meanipcconf1=0.0;
		double meanipcconf2=0.0;
		double meanipcconf3=0.0;
		double meanipcconf4=0.0;
		double meanipcconf5=0.0;
		double meanipcstable=0.0;
		double meaninstconf0=0.0;
		double meaninstconf1=0.0;
		double meaninstconf2=0.0;
		double meaninstconf3=0.0;
		double meaninstconf4=0.0;
		double meaninstconf5=0.0;
		int numberconf0=0;
		int numberconf1=0;
		int numberconf2=0;
		int numberconf3=0;
		int numberconf4=0;
		int numberconf5=0;
		int numberipcstable=0;
	

	std::set<std::string>::iterator it;
        //for each label
        for (it=_wdDescriptionSet.begin(); it!=_wdDescriptionSet.end(); ++it){

                std::pair <std::multimap<std::string ,nanos::WorkDescriptor::wdinfo >::iterator, std::multimap<std::string,nanos::WorkDescriptor::wdinfo>::iterator> ret;

		ret = _hashSysStatistics.equal_range(*it);

		//std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::iterator it3=ret.first;
		//if (it3!=NULL) myfile << getNumThreads()<<":"<< it3->second.slabel;
		
		int size= _hashSysStatistics.size();

		meanipcconf0=0.0;
		meanipcconf1=0.0;
		meanipcconf2=0.0;
		meanipcconf3=0.0;
		meanipcconf4=0.0;
		meanipcconf5=0.0;
		meanipcstable=0.0;
		meaninstconf0=0.0;
		meaninstconf1=0.0;
		meaninstconf2=0.0;
		meaninstconf3=0.0;
		meaninstconf4=0.0;
		meaninstconf5=0.0;
		numberconf0=0;
		numberconf1=0;
		numberconf2=0;
		numberconf3=0;
		numberconf4=0;
		numberconf5=0;
		numberipcstable=0;
	
		if(size!=0){

			bool cent=true;
			for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::iterator it2=ret.first; it2!=ret.second; ++it2){
				if(cent){ 
					fileipc << getNumThreads() << ":" <<  it2->second.wdlabel; cent=false;
					fileinst << getNumThreads() << ":" <<  it2->second.wdlabel; cent=false;
					//std::cout<<"task: "<<it2->second.wdlabel<< std::endl;
					}
	
				//fileipc << ":" << it2->second.ipc;
				//fileinst << ":" << it2->second.ins;
				
				if(it2->second.prefetcherconf==0 && it2->second.phasemode==0){
						meanipcconf0+=it2->second.ipc;
						meaninstconf0+=it2->second.ins;
						numberconf0+=1;
				}else if(it2->second.prefetcherconf==1 && it2->second.phasemode==0){
						meanipcconf1+=it2->second.ipc;
						meaninstconf1+=it2->second.ins;
						numberconf1+=1;
				}else if(it2->second.prefetcherconf==2 && it2->second.phasemode==0){
						meanipcconf2+=it2->second.ipc;
						meaninstconf2+=it2->second.ins;
						numberconf2+=1;
				}else if(it2->second.prefetcherconf==3 && it2->second.phasemode==0){
						meanipcconf3+=it2->second.ipc;
						meaninstconf3+=it2->second.ins;
						numberconf3+=1;
				}else if(it2->second.prefetcherconf==4 && it2->second.phasemode==0){
						meanipcconf4+=it2->second.ipc;
						meaninstconf4+=it2->second.ins;
						numberconf4+=1;
				}else if(it2->second.prefetcherconf==5 && it2->second.phasemode==0){
						meanipcconf5+=it2->second.ipc;
						meaninstconf5+=it2->second.ins;
						numberconf5+=1;
				}else if(it2->second.phasemode==1){
						meanipcstable+=it2->second.ipc;
						numberipcstable+=1;
				}	
				
			}

		/*std::cout<<"meanipcconf0: "<<meanipcconf0 <<" "<<" numberconf0: "<<numberconf0<< " division "<< float(meanipcconf0/float(numberconf0)) <<std::endl;
		std::cout<<"meanipcconf1: "<<meanipcconf1 <<" "<<" numberconf1: "<<numberconf1<< " division "<< float(meanipcconf1/float(numberconf1)) <<std::endl;
		std::cout<<"meanipcconf2: "<<meanipcconf2 <<" "<<" numberconf2: "<<numberconf2<< " division "<< float(meanipcconf2/float(numberconf2)) <<std::endl;
		std::cout<<"meanipcconf3: "<<meanipcconf3 <<" "<<" numberconf3: "<<numberconf3<< " division "<< float(meanipcconf3/float(numberconf3)) <<std::endl;
		std::cout<<"meanipcconf4: "<<meanipcconf4 <<" "<<" numberconf4: "<<numberconf4<< " division "<< float(meanipcconf4/float(numberconf4)) <<std::endl;
		std::cout<<"meanipcconf5: "<<meanipcconf5 <<" "<<" numberconf4: "<<numberconf5<< " division "<< float(meanipcconf5/float(numberconf5)) <<std::endl;*/

		meanipcconf0/=(float)numberconf0;
		if(numberconf1!=0){
			meanipcconf1/=(float)numberconf1;}
		else{
			meanipcconf1=meanipcconf0;
		}
		if(numberconf2!=0){
			meanipcconf2/=(float)numberconf2;}
		else{
			meanipcconf2=(meanipcconf0+meanipcconf1) /2.0;
		}
		if(numberconf3!=0){
			meanipcconf3/=(float)numberconf3;}
		else{
			meanipcconf3=(meanipcconf0+meanipcconf1+meanipcconf2) /3.0;
		}
		if(numberconf4!=0){
			meanipcconf4/=(float)numberconf4;}
		else{
			meanipcconf4=(meanipcconf0+meanipcconf1+meanipcconf2+meanipcconf3) /4.0;
		}
		if(numberconf5!=0){
			meanipcconf5/=(float)numberconf5;}
		else{
			meanipcconf5=(meanipcconf0+meanipcconf1+meanipcconf2+meanipcconf3+meanipcconf4) /5.0;
		}
		if(numberipcstable!=0){
			meanipcstable/=numberipcstable;}
		else{
			meanipcstable=meanipcconf5;
		}
		
		meaninstconf0/=(float)numberconf0;
		meaninstconf1/=(float)numberconf1;
		meaninstconf2/=(float)numberconf2;
		meaninstconf3/=(float)numberconf3;
		meaninstconf4/=(float)numberconf4;
		meaninstconf5/=(float)numberconf5;
		
//printf( "test %f\n", meanipcconf0);
//char * pointer;
//asprintf (&pointer, "%f",meanipcconf0);

		//fileipc <<":"<< *(new std::string (pointer)) <<":"<<(float)meanipcconf1<<":"<<(std::string)meanipcconf2<<":"<<(double)meanipcconf3<<":"<<meanipcconf4<<":"<<meanipcconf5<< std::endl;
		fileipc <<":"<< meanipcconf0 <<":"<<meanipcconf1<<":"<<meanipcconf2<<":"<<meanipcconf3<<":"<<meanipcconf4<<":"<<meanipcconf5<< std::endl;

		fileipc << numberconf0 <<":"<<numberconf1<<":"<<numberconf2<<":"<<numberconf3<<":"<<numberconf4<<":"<<numberconf5<< std::endl;
	
		fileipcstable <<":"<<meanipcstable<<std::endl;
		fileipcstable << numberipcstable << std::endl;

		fileinst <<":"<< meaninstconf0<<":"<<meaninstconf1<<":"<<meaninstconf2<<":"<<meaninstconf3<<":"<<meaninstconf4<<":"<<meaninstconf5<< std::endl;
		//fileipc  << std::endl;
		//fileinst << std::endl;
		}
	}
	//fileipc << getNumThreads()<<":ipc:"<<meantotalipc<< std::endl;
		
}

void System::writeDumpNumTimesConfs(){
	std::cout << "DumpingNumTimesConfs3"<< std::endl;
	std::ofstream filetimes;
	
	char threads[15]; sprintf(threads, "%d", _numThreads);
	char numps[15]; sprintf(numps, "%d", _numPS);
	char numpb[15]; sprintf(numpb, "%d", _numPB);
  char dot[1]; dot[0]='.'; char strtimes[61]; 
	strcpy (strtimes,threads); strcat (strtimes,dot); strcat (strtimes,numps); strcat (strtimes,dot); strcat (strtimes,numpb);
	strcat (strtimes,dot); //strcat (str,"statistics");

  if(_kbPowerBW>0) {
                char kpbw [15];//knobpowerbandwidth
                sprintf(kpbw, "%d", _kbPowerBW);
								strcat(strtimes,kpbw);strcat(strtimes,dot);}
	if(_numIns>0) {
                char numinstance [15];
                sprintf(numinstance, "%d", _numIns);
								strcat(strtimes,numinstance);strcat(strtimes,dot);}
	strcat(strtimes,"CONFS.txt");
  filetimes.open  ((char*)strtimes, std::ios::trunc );//app
	////////////////****
	int vnconfs [6];
	for (int i=0;i<6;i++){
		vnconfs[i]=0;
	}
	std::set<std::string>::iterator it; //for each label
        for (it=_wdDescriptionSet.begin(); it!=_wdDescriptionSet.end(); ++it){
                std::pair <std::multimap<std::string ,nanos::WorkDescriptor::wdinfo >::iterator, std::multimap<std::string,nanos::WorkDescriptor::wdinfo>::iterator> ret;
		ret = _hashSysStatistics.equal_range(*it);
		int size= _hashSysStatistics.size();
		if(size!=0){
			bool cent=true;
			for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::iterator it2=ret.first; it2!=ret.second; ++it2){
				if(cent){ 
					vnconfs[it2->second.prefetcherconf]+=1;
					}
			}
		}
	}
	for (int i=0;i<6;i++){
		filetimes << vnconfs[i] <<std::endl;
	}
}


void System::writeDumpOrderConfs(){
  std::cout << "DumpingOrdsConfs"<< std::endl;
  std::ofstream filetimes;

  char threads[15]; sprintf(threads, "%d", _numThreads);
  char numps[15]; sprintf(numps, "%d", _numPS);
  char numpb[15]; sprintf(numpb, "%d", _numPB);
  char dot[1]; dot[0]='.'; char strtimes[61];
  strcpy (strtimes,threads); strcat (strtimes,dot); strcat (strtimes,numps); strcat (strtimes,dot); strcat (strtimes,numpb);
  strcat (strtimes,dot); //strcat (str,"statistics");

  if(_kbPowerBW>0) {
                char kpbw [15];//knobpowerbandwidth
                sprintf(kpbw, "%d", _kbPowerBW);
                strcat(strtimes,kpbw);strcat(strtimes,dot);}
  if(_numIns>0) {
                char numinstance [15];
                sprintf(numinstance, "%d", _numIns);
                strcat(strtimes,numinstance);strcat(strtimes,dot);}
  strcat(strtimes,"ORDCONFS.txt");
  filetimes.open  ((char*)strtimes, std::ios::trunc );//app
  ////////////////****
  /*int vnconfs [6];
  for (int i=0;i<6;i++){
    vnconfs[i]=0;
  }*/
  std::set<std::string>::iterator it; //for each label
 for (it=_wdDescriptionSet.begin(); it!=_wdDescriptionSet.end(); ++it){
                std::pair <std::multimap<std::string ,nanos::WorkDescriptor::wdinfo >::iterator, std::multimap<std::string,nanos::WorkDescriptor::wdinfo>::iterator> ret;
    ret = _hashSysStatistics.equal_range(*it);
    int size= _hashSysStatistics.size();
    if(size!=0){
      bool cent=true;
      for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::iterator it2=ret.first; it2!=ret.second; ++it2){//for each task of a task type
        if(cent){

          filetimes << getNumThreads() << ":" <<  it2->second.wdlabel; cent=false;
          //vnconfs[it2->second.prefetcherconf]+=1;
          }
          if( it2->second.phasemode==1 ){
                filetimes <<":"<< it2->second.prefetcherconf;
          }
      }
      filetimes  << std::endl;
    }
  }
  //for (int i=0;i<6;i++){
  //  filetimes << vnconfs[i] <<std::endl;
  //}
}



///////////////////
void System::getStatisticsGroups(){
//std::cout << "getting statistics groups :\n";
  	std::set<std::string>::iterator it;
	//for each label
		int numtasks=0;
  	for (it=_wdDescriptionSet.begin(); it!=_wdDescriptionSet.end(); ++it){

    		std::pair <std::multimap<std::string ,nanos::WorkDescriptor::wdinfo >::iterator, std::multimap<std::string,nanos::WorkDescriptor::wdinfo>::iterator> ret;

    ret = _hashSysStatistics.equal_range(*it);
    //std::cout << *it << " =>";

		labelstats aux;//=new labelstats();
		
		aux.slabel= *it;
		aux.meanipc=0.0, aux.stdipc = 0.0 ;
	  aux.meanins=0.0, aux.stdins = 0.0;
	        aux.meandins=0.0, aux.stddins = 0.0;
        	aux.meanrtime=0.0, aux.stdrtime = 0.0;
         	aux.meanptime=0.0, aux.stdptime = 0.0;
         	aux.meanl1tcm=0.0, aux.stdl1tcm = 0.0;
        	aux.meanl2tcm=0.0, aux.stdl2tcm = 0.0;
         	aux.meanl3tcm=0.0, aux.stdl3tcm = 0.0;

		//for each task instance of that label
		
		long long numwd=0;
		//aux.minipc= 99999.0 ;
		//aux.maxipc= 0.0  ;
    		for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::iterator it2=ret.first; it2!=ret.second; ++it2){
			//if(aux.maxipc < it->second.ipc) aux.maxipc=it->second.ipc;
			//if(aux.minipc > it->second.ipc) aux.minipc=it->second.ipc;
			
			aux.meanipc+=it2->second.ipc;
			aux.meanins+=it2->second.ins;
			aux.meandins+=it2->second.dins;
			aux.meanrtime+=it2->second.rtime;
			aux.meanptime+=it2->second.ptime;
			aux.meanl1tcm+=it2->second.l1tcm;
			aux.meanl2tcm+=it2->second.l2tcm;
			aux.meanl3tcm+=it2->second.l3tcm;
			numwd++;
		}
								aux.meanipc   /= numwd;
                aux.meanins   /= numwd;
                aux.meandins  /= numwd;
                aux.meanrtime /= numwd;
                aux.meanptime /= numwd;
                aux.meanl1tcm /= numwd;
                aux.meanl2tcm /= numwd;
                aux.meanl3tcm /= numwd;
								meantotalipc+=aux.meanipc;	
								numtasks++;
		//for each task instance of that label
		for (std::multimap<std::string, nanos::WorkDescriptor::wdinfo>::iterator it2=ret.first; it2!=ret.second; ++it2){
			aux.stdipc= (it2->second.ipc - aux.meanipc)*(it2->second.ipc - aux.meanipc);
                	aux.stdins= (it2->second.ins - aux.meanins)*(it2->second.ins - aux.meanins);
                	aux.stddins= (it2->second.dins - aux.meandins)*(it2->second.dins - aux.meandins);
                	aux.stdrtime= (it2->second.rtime - aux.meanrtime)*(it2->second.rtime - aux.meanrtime);
                	aux.stdptime= (it2->second.ptime - aux.meanptime)*(it2->second.ptime - aux.meanptime);
                	aux.stdl1tcm= (it2->second.l1tcm - aux.meanl1tcm)*(it2->second.l1tcm - aux.meanl1tcm);
                	aux.stdl2tcm= (it2->second.l2tcm - aux.meanl2tcm)*(it2->second.l2tcm - aux.meanl2tcm);
                	aux.stdl2tcm= (it2->second.l3tcm - aux.meanl3tcm)*(it2->second.l3tcm - aux.meanl3tcm);
		}
		aux.stdipc  =sqrt(aux.stdipc/numwd);
	        aux.stdins  =sqrt(aux.stdins/numwd);
	        aux.stddins =sqrt(aux.stddins/numwd);
        	aux.stdrtime=sqrt(aux.stdrtime/numwd);
        	aux.stdptime=sqrt(aux.stdptime/numwd);
        	aux.stdl1tcm=sqrt(aux.stdl1tcm/numwd);
        	aux.stdl2tcm=sqrt(aux.stdl2tcm/numwd);
        	aux.stdl3tcm=sqrt(aux.stdl3tcm/numwd);
	
    		_sysStatisticsGroups.push_back(aux);
	
    		//std::cout << '\n';
  	}
	meantotalipc=meantotalipc/(float)(numtasks);
	//std::cout<<"meanipcs: "<< meantotalipc << std::endl;
}

void System::getStatistics(){
//	printStatistics();
//	writeStatistics();
//	printGroups();
	getStatisticsGroups();
////	writeStatisticsGroups();
	writeDump();
	writeDumpConfs();
	

writeDumpNumTimesConfs();	
  writeDumpOrderConfs();
}


//prints stdout merged task type statistics
void System::printStatistics (){
	std::cout << "Mean IPC: " << meanipc <<"\n";
	std::cout << "Mean Instructions: " << meanins <<"\n";
	std::cout << "Mean Real Time: " << meanrtime <<"\n";
	std::cout << "Mean Process Time: " << meanptime <<"\n";
	std::cout << "STD IPC: " << stdipc <<"\n";
	std::cout << "STD Instructions: " << stdins  <<"\n";
	std::cout << "STD Real Time: " << stdrtime <<"\n";
	std::cout << "STD process Time: " << stdptime <<"\n";
}
//writes to file merged task type statistics
void System::writeStatistics(){
	
	std::ofstream myfile;
	myfile.open ("ompssstatistics.txt", std::ios::app);

	myfile <<getNumThreads()<<":"<<getCreatedTasks()<<":"<<"0"<<":"<<meanipc<<":"<<stdipc<<":"<<meanins<<":"<<stdins<<":"<<meanrtime<<":"<<stdrtime<<":"<<meanptime<<":"<<stdptime<<":"<<meanl1tcm<<":"<<stdl1tcm<<":"<<meanl2tcm<<":"<<stdl2tcm<<":"<<meanl3tcm<<":"<<stdl3tcm<< "\n"; 

	myfile.close();	

	std::ofstream myfileh;
	myfileh.open ("ompssshuman", std::ios::app);
	myfileh <<"NThreads:"<<getNumThreads()<<" "<<"XIPC:"<<meanipc<<" "<<"XL1TCM:"<<meanl1tcm<<" "<<"XL2TCM:"<<meanl2tcm<<" "<<"XL3TCM:"<<meanl3tcm<<"\n";
	
	myfileh.close();
}

//writes to file for each task type its generalized statistics
void System::writeStatisticsGroups(){
	std::cout << "statistcisGroups:" << "\n";
	std::ofstream myfile;
	//std::stringstream ss;
	//ss << _numThreads;

	//std::string ompssstatisticshash= ss.str().c_str();//+"."+_numPS+"."+_numPB;
	//char exec[15];
	//sprintf(exec, "%d", argv[0]);
	
	char threads[15]; sprintf(threads, "%d", _numThreads);
	char numps[15]; sprintf(numps, "%d", _numPS);
	char numpb[15]; sprintf(numpb, "%d", _numPB);
  char dot[1]; dot[0]='.'; char str[61];
	//strcpy (str,exec);
	//strcat (str,dot);
	strcpy (str,threads); strcat (str,dot); strcat (str,numps); strcat (str,dot); strcat (str,numpb);
	strcat (str,dot); strcat (str,"statistics");

	myfile.open((char*)str, std::ios::app);
	
	for (std::list<labelstats>::iterator ci = _sysStatisticsGroups.begin(); ci != _sysStatisticsGroups.end(); ++ci){
	//	std::cout << "label: " << ci->slabel << "\n";
	//	myfile <<  ci->slabel <<":"<<ci->meanipc":"<<ci->stdipc<<":"<<ci->meanins<<":"<<ci->stdins<<":"<<ci->rtime <<":"<<ci-> <<":"<<   <<"\n" ;

	myfile <<getNumThreads()<<":"<<ci->slabel<<":"<<ci->meanipc<<":"<<ci->stdipc<<":"<<ci->meanins<<":"<<ci->stdins<<":"<<ci->meanrtime<<":"<<ci->stdrtime<<":"<<ci->meanptime<<":"<<ci->stdptime<<":"<<ci->meanl1tcm<<":"<<ci->stdl1tcm<<":"<<ci->meanl2tcm<<":"<<ci->stdl2tcm<<":"<<ci->meanl3tcm<<":"<<ci->stdl3tcm<< "\n"; 

	}
	myfile.close();
}

/*void System::writeStatisticsPlots(){
	std::ofstream labels, threads, values;
	labels.open("labels.txt", std::ios::app);
	threads.open("labels.txt", std::ios::app);
	values.open("labels.txt", std::ios::app);
	

	for (std::list<labelstats>::iterator ci = _sysStatisticsGroups.begin(); ci != _sysStatisticsGroups.end(); ++ci){
		labels << "'"<< ci->slabel <<"',";
						 
	}

}*/

#endif

std::vector<std::string> &System::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> System::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


// default system values go here 
System::System () :
      _atomicWDSeed( 1 ), _threadIdSeed( 0 ),
      _numPEs( INT_MAX ), _numPS(0), _numPB(0), _kbPowerBW(0) , _applyTrack(0), _numThreads( 0 ), _deviceStackSize( 0 ), _bindingStart (0), _bindingStride(1),  _bindThreads( true ), _profile( false ),
      _instrument( false ), _verboseMode( false ), _summary( false ), _executionMode( DEDICATED ), _initialMode( POOL ),
      _untieMaster( true ), _delayedStart( false ), _useYield( true ), _synchronizedStart( true ),
      _numSockets( 0 ), _coresPerSocket( 0 ), _numAvailSockets( 0 ), _enable_dlb( false ), _throttlePolicy ( NULL ),
      _schedStats(), _schedConf(), _defSchedule( "bf" ), _defThrottlePolicy( "hysteresis" ), 
      _defBarr( "centralized" ), _defInstr ( "empty_trace" ), _defDepsManager( "plain" ), _defArch( "smp" ),
      _initializedThreads ( 0 ), _targetThreads ( 0 ), _pausedThreads( 0 ),
      _pausedThreadsCond(), _unpausedThreadsCond(),
      _instrumentation ( NULL ), _defSchedulePolicy( NULL ), _dependenciesManager( NULL ),
      _pmInterface( NULL ), _useCaches( true ), _cachePolicy( System::DEFAULT ), _cacheMap()

#ifdef GPU_DEV
      , _pinnedMemoryCUDA( NEW CUDAPinnedMemoryManager() )
#endif
#ifdef NANOS_INSTRUMENTATION_ENABLED
      , _enableEvents(), _disableEvents(), _instrumentDefault("default"), _enable_cpuid_event( false )
#endif
      , _lockPoolSize(37), _lockPool( NULL ), _mainTeam (NULL)
{
   verbose0 ( "NANOS++ initializing... start" );

		
	 meantotalipc=0.0; //davidp
   // OS::init must be called here and not in System::start() as it can be too late
   // to locate the program arguments at that point

   OS::init();
   config();

	 std::cout << "START: " << _applyTrack  << std::endl;
	 std::cout << "value of numInstance: " << _numIns  << std::endl;
	 std::cout << "value of knob power bw: " << _kbPowerBW  << std::endl;
	 std::cout << "value of tracking: " << _applyTrack  << std::endl;

	
   OS::getProcessAffinity( &_cpu_set );

		if(_numPS==0) {	_numPS=16;}
		if(_numPB==0) { _numPB=4000;}
		if(_kbPowerBW==0) {_kbPowerBW=0;}
	 
		///if required, fill the track with a file////////////////////!!!
		if(_applyTrack==1){
				std::cout << "DumpingOrdsConfs _numThreads:"<< _numThreads <<std::endl;
				//std::ifstream filetimes;
				char threads[15]; sprintf(threads, "%d", _numPEs);
				char numps[15]; sprintf(numps, "%d", _numPS);
				char numpb[15]; sprintf(numpb, "%d", _numPB);
				char dot[1]; dot[0]='.'; char strtimes[61];
				strcpy (strtimes,threads); strcat (strtimes,dot); strcat (strtimes,numps); strcat (strtimes,dot); strcat (strtimes,numpb);
				strcat (strtimes,dot); //strcat (str,"statistics");

				if(_kbPowerBW>0) {
											char kpbw [15];//knobpowerbandwidth
											sprintf(kpbw, "%d", _kbPowerBW);
											strcat(strtimes,kpbw);strcat(strtimes,dot);}
				if(_numIns>0) {
											char numinstance [15];
											sprintf(numinstance, "%d", _numIns);
											strcat(strtimes,numinstance);strcat(strtimes,dot);}
				strcat(strtimes,"ORDCONFS.txt");
				
				std::ifstream filetimes (strtimes);
				//filetimes.open  ((char*)strtimes, std::ios::trunc );//app
				if (filetimes.is_open()){
					 std::string line;
					 std::string thislabel;
					 while ( getline (filetimes,line) ){
						//std::cout<< line <<std::endl		;
						//std::vector<std::string> x; //= std::split(line, ':');
						///std::stringstream ss(line);
						///std::istream_iterator<std::string> begin(ss);
						///std::istream_iterator<std::string> end;
						///std::vector<std::string> x(begin, end);
						///std::copy(x.begin(), x.end(), std::ostream_iterator<std::string>(std::cout, ":"));
						////std::vector<std::string> x = split(line, ':');
						x = split(line, ':');
						////std::vector<int> ints;
						////std::vector<std::string>::iterator itstr;
						////std::vector<int>::iterator itint;
						itstr=x.begin();
						itint=ints.begin();
						itstr++;
						thislabel=*itstr;
						//->std::cout << "label is: " << thislabel << std::endl;
						itstr++;
						
						/*for (; itstr<x.end();itstr++ ){
								std::cout	<< *itstr <<".";
						}
						std::cout<< std::endl;*/

		
					//->std::cout << "before copying without first two positions: "<< std::endl;
						for(  ; itstr<x.end(); itstr++  ){
							itint=ints.insert(itint	, atoi( (*itstr).c_str() ) );
							//ints.insert(itint	, 2 );
							//std::cout<<"Inserted an element: "<< atoi((*itstr).c_str())  <<std::endl;
							//itint++ ;
      				//cout << line << '\n';
    				}
						
					//->std::cout << "DATA COPIED: "<< std::endl;
					/*	for (itint=ints.begin(); itint<ints.end();itint++ ){
								std::cout	<< *itint <<".";
						}
						std::cout<< std::endl;*/
					}
						

					_hashSysTracking.insert( std::pair< std::string, std::vector<int> > (thislabel , ints));
					_hashSysCountTracking.insert( std::pair< std::string, int > (thislabel ,0 ));

    			filetimes.close();			
					std::cout << "endl transfering files to vector: "<< std::endl;
				}else{
						std::cout<<"unable to open: " << strtimes << std::endl;
				}
			//	_hasSysTracking  
		}
	

 
   int cpu_count = getCpuCount();

   std::vector<int> cpu_affinity;
   cpu_affinity.reserve( cpu_count );
   std::ostringstream oss_cpu_idx;
   oss_cpu_idx << "[";
   for ( int i=0; i<CPU_SETSIZE; i++ ) {
      if ( CPU_ISSET(i, &_cpu_set) ) {
         cpu_affinity.push_back(i);
         oss_cpu_idx << i << ", ";
      }
   }
   oss_cpu_idx << "]";
   
   verbose0("PID[" << getpid() << "]. CPU affinity " << oss_cpu_idx.str());
   
   // Ensure everything is properly configured
   if( getNumPEs() == INT_MAX && _numThreads == 0 )
      // If no parameter specified, use all available CPUs
      setNumPEs( cpu_count );
   
   if ( _numThreads == 0 )
      // No threads specified? Use as many as PEs
      _numThreads = _numPEs;
   else if ( getNumPEs() == INT_MAX ){
      // No number of PEs given? Use 1 thread per PE
      setNumPEs(  _numThreads );
   }

   // Set _bindings structure once we have the system mask and the binding info
   _bindings.reserve( cpu_count );
   for ( int i=0, collisions = 0; i < cpu_count; ) {

      // The cast over cpu_affinity is needed because std::vector::size() returns a size_t type
      int pos = (_bindingStart + _bindingStride*i + collisions) % (int)cpu_affinity.size();

      // 'pos' may be negative if either bindingStart or bindingStride were negative
      // this loop fixes that % operator is the remainder, not the modulo operation
      while ( pos < 0 ) pos+=cpu_affinity.size();

      if ( std::find( _bindings.begin(), _bindings.end(), cpu_affinity[pos] ) != _bindings.end() ) {
         collisions++;
         ensure( collisions != cpu_count, "Reached limit of collisions. We should never get here." );
         continue;
      }
      _bindings.push_back( cpu_affinity[pos] );
      i++;
   }

   CPU_ZERO( &_cpu_active_set );

   _lockPool = NEW Lock[_lockPoolSize];

   if ( !_delayedStart ) {
      start();
   }

////////////davidp
/*
	std::cout<<"PAPI ERROR before initialization: \n";
	int retval=0;
	if (not PAPI_is_initialized()) {

                if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT ){
                        std::cout<<"PAPI ERROR init: " << retval << "\n";
                }
                if((retval = PAPI_thread_init(pthread_self)) != PAPI_OK){
                        std::cout<<"PAPI ERROR thread init: " << retval << "\n";
                }

        }else{
                std::cout<<"PAPI ERROR already initialized: \n";
        }
*/
//\davidp
   verbose0 ( "NANOS++ initializing... end" );
}

struct LoadModule
{
   void operator() ( const char *module )
   {
      if ( module ) {
        verbose0( "loading " << module << " module" );
        sys.loadPlugin(module);
      }
   }
};

void System::loadModules ()
{
   verbose0 ( "Configuring module manager" );

   _pluginManager.init();

   verbose0 ( "Loading modules" );

   const OS::ModuleList & modules = OS::getRequestedModules();
   std::for_each(modules.begin(),modules.end(), LoadModule());

   // load host processor module
   if ( _hostFactory == NULL ) {
     verbose0( "loading Host support" );

     if ( !loadPlugin( "pe-"+getDefaultArch() ) )
       fatal0 ( "Couldn't load host support" );
   }
   ensure0( _hostFactory,"No default host factory" );

#ifdef GPU_DEV
   verbose0( "loading GPU support" );

   if ( !loadPlugin( "pe-gpu" ) )
      fatal0 ( "Couldn't load GPU support" );
#endif
   
#ifdef OpenCL_DEV
   verbose0( "loading OpenCL support" );

   if ( !loadPlugin( "pe-opencl" ) )
     fatal0 ( "Couldn't load OpenCL support" );
#endif
   
   // load default schedule plugin
   verbose0( "loading " << getDefaultSchedule() << " scheduling policy support" );

   if ( !loadPlugin( "sched-"+getDefaultSchedule() ) )
      fatal0 ( "Couldn't load main scheduling policy" );

   ensure0( _defSchedulePolicy,"No default system scheduling factory" );

   verbose0( "loading " << getDefaultThrottlePolicy() << " throttle policy" );

   if ( !loadPlugin( "throttle-"+getDefaultThrottlePolicy() ) )
      fatal0( "Could not load main cutoff policy" );

   ensure0( _throttlePolicy, "No default throttle policy" );

   verbose0( "loading " << getDefaultBarrier() << " barrier algorithm" );

   if ( !loadPlugin( "barrier-"+getDefaultBarrier() ) )
      fatal0( "Could not load main barrier algorithm" );

   if ( !loadPlugin( "instrumentation-"+getDefaultInstrumentation() ) )
      fatal0( "Could not load " + getDefaultInstrumentation() + " instrumentation" );

   ensure0( _defBarrFactory,"No default system barrier factory" );
   
   // load default dependencies plugin
   verbose0( "loading " << getDefaultDependenciesManager() << " dependencies manager support" );

   if ( !loadPlugin( "deps-"+getDefaultDependenciesManager() ) )
      fatal0 ( "Couldn't load main dependencies manager" );

   ensure0( _dependenciesManager,"No default dependencies manager" );

}

void System::unloadModules ()
{   
   delete _throttlePolicy;
   
   delete _defSchedulePolicy;
   
   //! \todo (#613): delete GPU plugin?
}

// Config Functor
struct ExecInit
{
   std::set<void *> _initialized;

   ExecInit() : _initialized() {}

   void operator() ( const nanos_init_desc_t & init )
   {
      if ( _initialized.find( (void *)init.func ) == _initialized.end() ) {
         init.func( init.data );
         _initialized.insert( ( void * ) init.func );
      }
   }
};

void System::config ()
{
   Config cfg;

   const OS::InitList & externalInits = OS::getInitializationFunctions();
   std::for_each(externalInits.begin(),externalInits.end(), ExecInit());
   
   if ( !_pmInterface ) {
      // bare bone run
      _pmInterface = NEW PMInterface();
   }

   //! Declare all configuration core's flags
   verbose0( "Preparing library configuration" );

   cfg.setOptionsSection( "Core", "Core options of the core of Nanos++ runtime" );

	 /////////////////////////////////////////
	 ////////////davidp test times prefetcher
	 cfg.registerConfigOption( "num-ps", NEW Config::IntegerVar( _numPS ),
                             "Number of times a prefetcher configuration is monitored" );
   cfg.registerArgOption( "num-ps", "num-ps" );
   cfg.registerEnvOption( "num-ps", "NX_NPS" );
	 //////////////davidp best times prefetcher
   cfg.registerConfigOption( "num-pb", NEW Config::IntegerVar( _numPB ),
                             "Number of times the best prefetcher configuration is used" );
   cfg.registerArgOption( "num-pb", "num-pb" );
   cfg.registerEnvOption( "num-pb", "NX_NPB" );
	 /////////////////////////////////

   cfg.registerConfigOption( "num_pes", NEW Config::UintVar( _numPEs ),
                             "Defines the number of processing elements" );
   cfg.registerArgOption( "num_pes", "pes" );
   cfg.registerEnvOption( "num_pes", "NX_PES" );
	///////////////////////////////////////
			  cfg.registerConfigOption( "num_ins", NEW Config::UintVar( _numIns ),
                             "Defines what number is an instance" );
        cfg.registerArgOption( "num_ins", "ins" );
        cfg.registerEnvOption( "num_ins", "NX_INS" );
	/////////////

		cfg.registerConfigOption( "apply_track", NEW Config::UintVar ( _applyTrack ),
                             "enables prfetcher conf tracking" );
   cfg.registerArgOption( "apply_track", "apt" );
   cfg.registerEnvOption( "apply_track", "NX_APPLYTRACK" );


				cfg.registerConfigOption( "knob_powerbw", NEW Config::UintVar( _kbPowerBW ),
																"power knob 1");
                             //"Knob to tune sensibility when finding most optimal configuration in terms of consumed power" );
        cfg.registerArgOption( "knob_powerbw", "kbpowerbw" );
        cfg.registerEnvOption( "knob_powerbw", "NX_KBBW" );

/*	 cfg.registerConfigOption( "apply-track", NEW Config::FlagOption( _applyTrack, false ),
                             "enables prfetcher conf tracking" );
   cfg.registerArgOption( "apply-track", "apt" );*/

	
	/////////////////
   cfg.registerConfigOption( "num_threads", NEW Config::PositiveVar( _numThreads ),
                             "Defines the number of threads. Note that OMP_NUM_THREADS is an alias to this." );
   cfg.registerArgOption( "num_threads", "threads" );
   cfg.registerEnvOption( "num_threads", "NX_THREADS" );
   
   cfg.registerConfigOption( "cores-per-socket", NEW Config::PositiveVar( _coresPerSocket ),
                             "Number of cores per socket." );
   cfg.registerArgOption( "cores-per-socket", "cores-per-socket" );
   
   cfg.registerConfigOption( "num-sockets", NEW Config::PositiveVar( _numSockets ),
                             "Number of sockets available." );
   cfg.registerArgOption( "num-sockets", "num-sockets" );
   
   cfg.registerConfigOption( "hwloc-topology", NEW Config::StringVar( _topologyPath ),
                             "Overrides hwloc's topology discovery and uses the one provided by an XML file." );
   cfg.registerArgOption( "hwloc-topology", "hwloc-topology" );
   cfg.registerEnvOption( "hwloc-topology", "NX_HWLOC_TOPOLOGY_PATH" );
   

   cfg.registerConfigOption( "stack-size", NEW Config::PositiveVar( _deviceStackSize ),
                             "Defines the default stack size for all devices" );
   cfg.registerArgOption( "stack-size", "stack-size" );
   cfg.registerEnvOption( "stack-size", "NX_STACK_SIZE" );

   cfg.registerConfigOption( "binding-start", NEW Config::IntegerVar ( _bindingStart ),
                             "Set initial cpu id for binding (binding required)" );
   cfg.registerArgOption( "binding-start", "binding-start" );
   cfg.registerEnvOption( "binding-start", "NX_BINDING_START" );

   cfg.registerConfigOption( "binding-stride", NEW Config::IntegerVar ( _bindingStride ),
                             "Set binding stride (binding required)" );
   cfg.registerArgOption( "binding-stride", "binding-stride" );
   cfg.registerEnvOption( "binding-stride", "NX_BINDING_STRIDE" );

   cfg.registerConfigOption( "no-binding", NEW Config::FlagOption( _bindThreads, false ),
                             "Disables thread binding" );
   cfg.registerArgOption( "no-binding", "disable-binding" );

   cfg.registerConfigOption( "no-yield", NEW Config::FlagOption( _useYield, false ),
                             "Do not yield on idle and condition waits" );
   cfg.registerArgOption( "no-yield", "disable-yield" );

   cfg.registerConfigOption( "verbose", NEW Config::FlagOption( _verboseMode ),
                             "Activates verbose mode" );
   cfg.registerArgOption( "verbose", "verbose" );

   cfg.registerConfigOption( "summary", NEW Config::FlagOption( _summary ),
                             "Activates summary mode" );
   cfg.registerArgOption( "summary", "summary" );

//! \bug implement execution modes (#146) */
#if 0
   cfg::MapVar<ExecutionMode> map( _executionMode );
   map.addOption( "dedicated", DEDICATED).addOption( "shared", SHARED );
   cfg.registerConfigOption ( "exec_mode", &map, "Execution mode" );
   cfg.registerArgOption ( "exec_mode", "mode" );
#endif

   registerPluginOption( "schedule", "sched", _defSchedule,
                         "Defines the scheduling policy", cfg );
   cfg.registerArgOption( "schedule", "schedule" );
   cfg.registerEnvOption( "schedule", "NX_SCHEDULE" );

   registerPluginOption( "throttle", "throttle", _defThrottlePolicy,
                         "Defines the throttle policy", cfg );
   cfg.registerArgOption( "throttle", "throttle" );
   cfg.registerEnvOption( "throttle", "NX_THROTTLE" );

   cfg.registerConfigOption( "barrier", NEW Config::StringVar ( _defBarr ),
                             "Defines barrier algorithm" );
   cfg.registerArgOption( "barrier", "barrier" );
   cfg.registerEnvOption( "barrier", "NX_BARRIER" );

   registerPluginOption( "instrumentation", "instrumentation", _defInstr,
                         "Defines instrumentation format", cfg );
   cfg.registerArgOption( "instrumentation", "instrumentation" );
   cfg.registerEnvOption( "instrumentation", "NX_INSTRUMENTATION" );

   cfg.registerConfigOption( "no-sync-start", NEW Config::FlagOption( _synchronizedStart, false),
                             "Disables synchronized start" );
   cfg.registerArgOption( "no-sync-start", "disable-synchronized-start" );

   cfg.registerConfigOption( "architecture", NEW Config::StringVar ( _defArch ),
                             "Defines the architecture to use (smp by default)" );
   cfg.registerArgOption( "architecture", "architecture" );
   cfg.registerEnvOption( "architecture", "NX_ARCHITECTURE" );

   cfg.registerConfigOption( "no-caches", NEW Config::FlagOption( _useCaches, false ), "Disables the use of caches" );
   cfg.registerArgOption( "no-caches", "disable-caches" );

   CachePolicyConfig *cachePolicyCfg = NEW CachePolicyConfig ( _cachePolicy );
   cachePolicyCfg->addOption("wt", System::WRITE_THROUGH );
   cachePolicyCfg->addOption("wb", System::WRITE_BACK );
   cachePolicyCfg->addOption( "nocache", System::NONE );

   cfg.registerConfigOption( "cache-policy", cachePolicyCfg,
                             "Defines the general cache policy to use: write-through / write-back. Can be overwritten for specific architectures" );
   cfg.registerArgOption( "cache-policy", "cache-policy" );
   cfg.registerEnvOption( "cache-policy", "NX_CACHE_POLICY" );
   
   registerPluginOption( "deps", "deps", _defDepsManager,
                         "Defines the dependencies plugin", cfg );
   cfg.registerArgOption( "deps", "deps" );
   cfg.registerEnvOption( "deps", "NX_DEPS" );
   

#ifdef NANOS_INSTRUMENTATION_ENABLED
   cfg.registerConfigOption( "instrument-default", NEW Config::StringVar ( _instrumentDefault ),
                             "Set instrumentation event list default (none, all)" );
   cfg.registerArgOption( "instrument-default", "instrument-default" );

   cfg.registerConfigOption( "instrument-enable", NEW Config::StringVarList ( _enableEvents ),
                             "Add events to instrumentation event list" );
   cfg.registerArgOption( "instrument-enable", "instrument-enable" );

   cfg.registerConfigOption( "instrument-disable", NEW Config::StringVarList ( _disableEvents ),
                             "Remove events to instrumentation event list" );
   cfg.registerArgOption( "instrument-disable", "instrument-disable" );

   cfg.registerConfigOption( "instrument-cpuid", NEW Config::FlagOption ( _enable_cpuid_event ),
                             "Add cpuid event when binding is disabled (expensive)" );
   cfg.registerArgOption( "instrument-cpuid", "instrument-cpuid" );
#endif

   cfg.registerConfigOption( "enable-dlb", NEW Config::FlagOption ( _enable_dlb ),
                              "Tune Nanos Runtime to be used with Dynamic Load Balancing library)" );
   cfg.registerArgOption( "enable-dlb", "enable-dlb" );

   _schedConf.config( cfg );
   _pmInterface->config( cfg );

   verbose0 ( "Reading Configuration" );

   cfg.init();
}

PE * System::createPE ( std::string pe_type, int pid, int uid )
{
   //! \todo lookup table for PE factories, in the mean time assume only one factory
   return _hostFactory( pid, uid );
}

void System::start ()
{
   if ( !_useCaches ) _cachePolicy = System::NONE;
   
   //! Load hwloc first, in order to make it available for modules
   if ( isHwlocAvailable() )
      loadHwloc();

   // loadNUMAInfo needs _targetThreads when hwloc is not available.
   // Note that it is not its final value!
   _targetThreads = _numThreads;
   
   // Load & check NUMA config
   loadNUMAInfo();

   // Modules can be loaded now
   loadModules();

   // Increase targetThreads, ask the architecture plugins
   for ( ArchitecturePlugins::const_iterator it = _archs.begin();
        it != _archs.end(); ++it )
   {
      _targetThreads += (*it)->getNumThreads();
   }

   // Instrumentation startup
   NANOS_INSTRUMENT ( sys.getInstrumentation()->filterEvents( _instrumentDefault, _enableEvents, _disableEvents ) );
   NANOS_INSTRUMENT ( sys.getInstrumentation()->initialize() );

   verbose0 ( "Starting runtime" );

   _pmInterface->start();

   int numPes = getNumPEs();

   _pes.reserve ( numPes );

   PE *pe = createPE ( _defArch, getBindingId( 0 ), 0 );
   pe->setNUMANode( getNodeOfPE( pe->getId() ) );
   _pes.push_back ( pe );
   _workers.push_back( &pe->associateThisThread ( getUntieMaster() ) );
   CPU_SET( getBindingId( 0 ), &_cpu_active_set );

   WD &mainWD = *myThread->getCurrentWD();
   (void) mainWD.getDirectory(true);
   
   if ( _pmInterface->getInternalDataSize() > 0 ) {
      char *data = NEW char[_pmInterface->getInternalDataSize()];
      _pmInterface->initInternalData( data );
      mainWD.setInternalData( data );
   }

   _pmInterface->setupWD( mainWD );

   /* Renaming currend thread as Master */
   myThread->rename("Master");

   NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_STARTUP) );
   
   // For each plugin, notify it's the way to reserve PEs if they are required
   for ( ArchitecturePlugins::const_iterator it = _archs.begin();
        it != _archs.end(); ++it )
   {
      (*it)->createBindingList();
   }   
   // Right now, _bindings should only store SMP PEs ids

   // Create PEs
   int p;
   for ( p = 1; p < numPes ; p++ ) {
      pe = createPE ( "smp", getBindingId( p ), p );
      pe->setNUMANode( getNodeOfPE( pe->getId() ) );
      _pes.push_back ( pe );

      CPU_SET( getBindingId( p ), &_cpu_active_set );
   }
   
   // Create threads
   for ( int ths = 1; ths < _numThreads; ths++ ) {
      pe = _pes[ ths % numPes ];
      _workers.push_back( &pe->startWorker() );
   }
   
   // For each plugin create PEs and workers
   //! \bug  FIXME (#855)
   for ( ArchitecturePlugins::const_iterator it = _archs.begin();
        it != _archs.end(); ++it )
   {
      for ( unsigned archPE = 0; archPE < (*it)->getNumPEs(); ++archPE )
      {
         PE * processor = (*it)->createPE( archPE, p );
         fatal_cond0( processor == NULL, "ArchPlugin::createPE returned NULL" );
         _pes.push_back( processor );
         _workers.push_back( &processor->startWorker() );
         CPU_SET( processor->getId(), &_cpu_active_set );
         ++p;
      }
   }

   // Set up internal data for each worker
   for ( ThreadList::const_iterator it = _workers.begin(); it != _workers.end(); it++ ) {

      WD & threadWD = (*it)->getThreadWD();
      if ( _pmInterface->getInternalDataSize() > 0 ) {
         char *data = NEW char[_pmInterface->getInternalDataSize()];
         _pmInterface->initInternalData( data );
         threadWD.setInternalData( data );
      }
      _pmInterface->setupWD( threadWD );
   }
      
#ifdef SPU_DEV
   PE *spu = NEW nanos::ext::SPUProcessor(100, (nanos::ext::SMPProcessor &) *_pes[0]);
   spu->startWorker();
#endif


   if ( getSynchronizedStart() ) {
      // Wait for the rest of the gang to be ready, but do not yet notify master thread is ready
      while (_initializedThreads.value() < ( _targetThreads - 1 ) ) {}
   }

   // FIXME (855): do this before thread creation, after PE creation
   completeNUMAInfo();

   switch ( getInitialMode() )
   {
      case POOL:
         verbose0("Pool model enabled (OmpSs)");
         _mainTeam = createTeam( _workers.size() );
         break;
      case ONE_THREAD:
         verbose0("One-thread model enabled (OpenMP)");
         _mainTeam = createTeam(1);
         break;
      default:
         fatal("Unknown initial mode!");
         break;
   }

   NANOS_INSTRUMENT ( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
   NANOS_INSTRUMENT ( static nanos_event_key_t num_threads_key = ID->getEventKey("set-num-threads"); )
   NANOS_INSTRUMENT ( nanos_event_value_t team_size =  (nanos_event_value_t) myThread->getTeam()->size(); )
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raisePointEvents(1, &num_threads_key, &team_size); )
   
   /* Master thread is ready now */
   if ( getSynchronizedStart() )
     threadReady();

   // Paused threads: set the condition checker 
   _pausedThreadsCond.setConditionChecker( EqualConditionChecker<unsigned int >( &_pausedThreads.override(), _workers.size() ) );
   _unpausedThreadsCond.setConditionChecker( EqualConditionChecker<unsigned int >( &_pausedThreads.override(), 0 ) );

   // All initialization is ready, call postInit hooks
   const OS::InitList & externalInits = OS::getPostInitializationFunctions();
   std::for_each(externalInits.begin(),externalInits.end(), ExecInit());

   NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent (NANOS_RUNNING) );
   
   // List unrecognised arguments
   std::string unrecog = Config::getOrphanOptions();
   if ( !unrecog.empty() )
      warning( "Unrecognised arguments: " << unrecog );
   Config::deleteOrphanOptions();
      
   // hwloc can be now unloaded
   if ( isHwlocAvailable() )
      unloadHwloc();

   if ( _summary )
      environmentSummary();
/////////////////////////////////////////////////////
#ifdef PAPI
//it has been execute before the associatethisthread method from processingelement
/*std::cout<<"PAPI ERROR before initialization (system) : \n";
        int retval=0;
        if (not PAPI_is_initialized()) {

                if((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT ){
                        std::cout<<"PAPI ERROR init (system): " << retval << "\n";
                }
                if((retval = PAPI_thread_init(pthread_self)) != PAPI_OK){
                        std::cout<<"PAPI ERROR thread init (system): " << retval << "\n";
                }

        }else{
                std::cout<<"PAPI ERROR already initialized (system): \n";
        }

std::cout<<"numthreads :  " << _numThreads << " num pes: " << numPes <<"\n";
*/
#endif

}

System::~System ()
{
   if ( !_delayedStart ) finish();
}

void System::finish ()
{
   /* Instrumentation: First removing RUNNING state from top of the state statck */
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseOpenStateEvent(NANOS_SHUTDOWN) );

   verbose ( "NANOS++ statistics");
   verbose ( std::dec << (unsigned int) getCreatedTasks() << " tasks has been executed" );

   verbose ( "NANOS++ shutting down.... init" );
   verbose ( "Wait for main workgroup to complete" );
   myThread->getCurrentWD()->waitCompletion( true );

   // we need to switch to the main thread here to finish
   // the execution correctly
   getMyThreadSafe()->getCurrentWD()->tied().tieTo(*_workers[0]);
   Scheduler::switchToThread(_workers[0]);
   
   ensure(getMyThreadSafe()->getId() == 0, "Main thread not finishing the application!");

   verbose ( "Joining threads... phase 1" );
   // signal stop PEs

   for ( unsigned p = 0; p < _pes.size() ; p++ ) {
      _pes[p]->stopAll();
   }

   verbose ( "Joining threads... phase 2" );

   // shutdown instrumentation
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raiseCloseStateEvent() );
   NANOS_INSTRUMENT ( sys.getInstrumentation()->finalize() );

   ensure( _schedStats._readyTasks == 0, "Ready task counter has an invalid value!");

   _pmInterface->finish();
   delete _pmInterface;

   /* System mem free */

   delete[] _lockPool;

   /* deleting master WD */
   //getMyThreadSafe()->getCurrentWD()->~WorkDescriptor();
   delete (WorkDescriptor *) (getMyThreadSafe()->getCurrentWD());

   for ( Slicers::const_iterator it = _slicers.begin(); it !=   _slicers.end(); it++ ) {
      delete (Slicer *)  it->second;
   }

   for ( WorkSharings::const_iterator it = _worksharings.begin(); it !=   _worksharings.end(); it++ ) {
      delete (WorkSharing *)  it->second;
   }
   
   /* deleting thread team */
   ThreadTeam* team = getMyThreadSafe()->getTeam();

   if ( team->getScheduleData() != NULL ) team->getScheduleData()->printStats();

   /* For every thread in the team */
   while ( team->size() > 0 ) {
      BaseThread* pThread = team->popThread();
      pThread->leaveTeam();
   }
   delete team;

   // join
   for ( unsigned p = 1; p < _pes.size() ; p++ ) {
      delete _pes[p];
   }
   
   /* unload modules */
   unloadModules();

   delete _dependenciesManager;

   // Deleting last processing element
   delete _pes[0];

   if ( allocator != NULL ) free (allocator);

   verbose ( "NANOS++ shutting down.... end" );

	_summary=true;
   if ( _summary )
      executionSummary();
}

/*! \brief Creates a new WD
 *
 *  This function creates a new WD, allocating memory space for device ptrs and
 *  data when necessary. 
 *
 *  \param [in,out] uwd is the related addr for WD if this parameter is null the
 *                  system will allocate space in memory for the new WD
 *  \param [in] num_devices is the number of related devices
 *  \param [in] devices is a vector of device descriptors 
 *  \param [in] data_size is the size of the related data
 *  \param [in,out] data is the related data (allocated if needed)
 *  \param [in] uwg work group to relate with
 *  \param [in] props new WD properties
 *  \param [in] num_copies is the number of copy objects of the WD
 *  \param [in] copies is vector of copy objects of the WD
 *  \param [in] num_dimensions is the number of dimension objects associated to the copies
 *  \param [in] dimensions is vector of dimension objects
 *
 *  When it does a full allocation the layout is the following:
 *  <pre>
 *  +---------------+
 *  |     WD        |
 *  +---------------+
 *  |    data       |
 *  +---------------+
 *  |  dev_ptr[0]   |
 *  +---------------+
 *  |     ....      |
 *  +---------------+
 *  |  dev_ptr[N]   |
 *  +---------------+
 *  |     DD0       |
 *  +---------------+
 *  |     ....      |
 *  +---------------+
 *  |     DDN       |
 *  +---------------+
 *  |    copy0      |
 *  +---------------+
 *  |     ....      |
 *  +---------------+
 *  |    copyM      |
 *  +---------------+
 *  |     dim0      |
 *  +---------------+
 *  |     ....      |
 *  +---------------+
 *  |     dimM      |
 *  +---------------+
 *  |   PM Data     |
 *  +---------------+
 *  </pre>
 */
void System::createWD ( WD **uwd, size_t num_devices, nanos_device_t *devices, size_t data_size, size_t data_align,
                        void **data, WG *uwg, nanos_wd_props_t *props, nanos_wd_dyn_props_t *dyn_props,
                        size_t num_copies, nanos_copy_data_t **copies, size_t num_dimensions,
                        nanos_region_dimension_internal_t **dimensions, nanos_translate_args_t translate_args,
                        const char *description )
{
   ensure(num_devices > 0,"WorkDescriptor has no devices");

   unsigned int i;
   char *chunk = 0;

   size_t size_CopyData;
   size_t size_Data, offset_Data, size_DPtrs, offset_DPtrs, size_Copies, offset_Copies, size_Dimensions, offset_Dimensions, offset_PMD;
   size_t offset_DESC, size_DESC;
   char *desc;
   size_t total_size;

   // WD doesn't need to compute offset, it will always be the chunk allocated address

   // Computing Data info
   size_Data = (data != NULL && *data == NULL)? data_size:0;
   if ( *uwd == NULL ) offset_Data = NANOS_ALIGNED_MEMORY_OFFSET(0, sizeof(WD), data_align );
   else offset_Data = 0; // if there are no wd allocated, it will always be the chunk allocated address

   // Computing Data Device pointers and Data Devicesinfo
   size_DPtrs    = sizeof(DD *) * num_devices;
   offset_DPtrs  = NANOS_ALIGNED_MEMORY_OFFSET(offset_Data, size_Data, __alignof__( DD*) );

   // Computing Copies info
   if ( num_copies != 0 ) {
      size_CopyData = sizeof(CopyData);
      size_Copies   = size_CopyData * num_copies;
      offset_Copies = NANOS_ALIGNED_MEMORY_OFFSET(offset_DPtrs, size_DPtrs, __alignof__(nanos_copy_data_t) );
      // There must be at least 1 dimension entry
      size_Dimensions = num_dimensions * sizeof(nanos_region_dimension_internal_t);
      offset_Dimensions = NANOS_ALIGNED_MEMORY_OFFSET(offset_Copies, size_Copies, __alignof__(nanos_region_dimension_internal_t) );
   } else {
      size_Copies = 0;
      // No dimensions
      size_Dimensions = 0;
      offset_Copies = offset_Dimensions = NANOS_ALIGNED_MEMORY_OFFSET(offset_DPtrs, size_DPtrs, 1);
   }

   // Computing description char * + description
   if ( description == NULL ) {
      offset_DESC = offset_Dimensions;
      size_DESC = size_Dimensions;
   } else {
      offset_DESC = NANOS_ALIGNED_MEMORY_OFFSET(offset_Dimensions, size_Dimensions, __alignof__ (void*) );
      size_DESC = (strlen(description)+1) * sizeof(char);
   }

   // Computing Internal Data info and total size
   static size_t size_PMD   = _pmInterface->getInternalDataSize();
   if ( size_PMD != 0 ) {
      static size_t align_PMD = _pmInterface->getInternalDataAlignment();
      offset_PMD = NANOS_ALIGNED_MEMORY_OFFSET(offset_DESC, size_DESC, align_PMD);
      total_size = NANOS_ALIGNED_MEMORY_OFFSET(offset_PMD,size_PMD,1);
   } else {
      offset_PMD = 0; // needed for a gcc warning
      total_size = NANOS_ALIGNED_MEMORY_OFFSET(offset_DESC, size_DESC, 1);
   }

   chunk = NEW char[total_size];
   if (props->clear_chunk)
       memset(chunk, 0, sizeof(char) * total_size);

   // allocating WD and DATA
   if ( *uwd == NULL ) *uwd = (WD *) chunk;
   if ( data != NULL && *data == NULL ) *data = (chunk + offset_Data);

   // allocating Device Data
   DD **dev_ptrs = ( DD ** ) (chunk + offset_DPtrs);
   for ( i = 0 ; i < num_devices ; i ++ ) dev_ptrs[i] = ( DD* ) devices[i].factory( devices[i].arg );

   ensure ((num_copies==0 && copies==NULL && num_dimensions==0 && dimensions==NULL) || (num_copies!=0 && copies!=NULL && num_dimensions!=0 && dimensions!=NULL ), "Number of copies and copy data conflict" );
   

   // allocating copy-ins/copy-outs
   if ( copies != NULL && *copies == NULL ) {
      *copies = ( CopyData * ) (chunk + offset_Copies);
      *dimensions = ( nanos_region_dimension_internal_t * ) ( chunk + offset_Dimensions );
   }

   // Copying description string
   if ( description == NULL ) desc = NULL;
   else {
      desc = (chunk + offset_DESC);
      strncpy ( desc, description, size_DESC);
//      desc[strlen(description)]='\0';
   }

   WD * wd =  new (*uwd) WD( num_devices, dev_ptrs, data_size, data_align, data != NULL ? *data : NULL,
                             num_copies, (copies != NULL)? *copies : NULL, translate_args, desc );
   // Set WD's socket
   wd->setSocket( getCurrentSocket() );
   
   // Set total size
   wd->setTotalSize(total_size );
   
   if ( getCurrentSocket() >= sys.getNumSockets() )
      throw NANOS_INVALID_PARAM;

   // All the implementations for a given task will have the same ID
   wd->setVersionGroupId( ( unsigned long ) devices );

   // initializing internal data
   if ( size_PMD > 0) {
      _pmInterface->initInternalData( chunk + offset_PMD );
      wd->setInternalData( chunk + offset_PMD );
   }

   // add to workgroup
   if ( uwg != NULL ) {
      WG * wg = ( WG * )uwg;
      wg->addWork( *wd );
   }

   // set properties
   if ( props != NULL ) {
      if ( props->tied ) wd->tied();
      wd->setPriority( dyn_props->priority );
      wd->setFinal ( dyn_props->flags.is_final );
   }
   if ( dyn_props && dyn_props->tie_to ) wd->tieTo( *( BaseThread * )dyn_props->tie_to );
   
   /* DLB */
   // In case the master have been busy crating tasks 
   // every 10 tasks created I'll check available cpus
   if(_atomicWDSeed.value()%10==0)dlb_updateAvailableCpus();
}

/*! \brief Creates a new Sliced WD
 *
 *  \param [in,out] uwd is the related addr for WD if this parameter is null the
 *                  system will allocate space in memory for the new WD
 *  \param [in] num_devices is the number of related devices
 *  \param [in] devices is a vector of device descriptors 
 *  \param [in] outline_data_size is the size of the related data
 *  \param [in,out] outline_data is the related data (allocated if needed)
 *  \param [in] uwg work group to relate with
 *  \param [in] slicer is the related slicer which contains all the methods to manage this WD
 *  \param [in,out] data used as the slicer data (allocated if needed)
 *  \param [in] props new WD properties
 *
 *  \return void
 *
 *  \par Description:
 * 
 *  This function creates a new Sliced WD, allocating memory space for device ptrs and
 *  data when necessary. Also allocates Slicer Data object which is related with the WD.
 *
 *  When it does a full allocation the layout is the following:
 *  <pre>
 *  +---------------+
 *  |   slicedWD    |
 *  +---------------+
 *  |    data       |
 *  +---------------+
 *  |  dev_ptr[0]   |
 *  +---------------+
 *  |     ....      |
 *  +---------------+
 *  |  dev_ptr[N]   |
 *  +---------------+
 *  |     DD0       |
 *  +---------------+
 *  |     ....      |
 *  +---------------+
 *  |     DDN       |
 *  +---------------+
 *  |    copy0      |
 *  +---------------+
 *  |     ....      |
 *  +---------------+
 *  |    copyM      |
 *  +---------------+
 *  |     dim0      |
 *  +---------------+
 *  |     ....      |
 *  +---------------+
 *  |     dimM      |
 *  +---------------+
 *  |   PM Data     |
 *  +---------------+
 *  </pre>
 *
 * \sa createWD, duplicateWD, duplicateSlicedWD
 */
void System::createSlicedWD ( WD **uwd, size_t num_devices, nanos_device_t *devices, size_t outline_data_size,
                        int outline_data_align, void **outline_data, WG *uwg, Slicer *slicer, nanos_wd_props_t *props,
                        nanos_wd_dyn_props_t *dyn_props, size_t num_copies, nanos_copy_data_t **copies, size_t num_dimensions,
                        nanos_region_dimension_internal_t **dimensions, const char *description )
{
   ensure(num_devices > 0,"WorkDescriptor has no devices");

   unsigned int i;
   char *chunk = 0;

   size_t size_CopyData;
   size_t size_Data, offset_Data, size_DPtrs, offset_DPtrs;
   size_t size_Copies, offset_Copies, size_Dimensions, offset_Dimensions, offset_PMD;
   size_t offset_DESC, size_DESC;
   char *desc;
   size_t total_size;

   // WD doesn't need to compute offset, it will always be the chunk allocated address

   // Computing Data info
   size_Data = (outline_data != NULL && *outline_data == NULL)? outline_data_size:0;
   if ( *uwd == NULL ) offset_Data = NANOS_ALIGNED_MEMORY_OFFSET(0, sizeof(SlicedWD), outline_data_align );
   else offset_Data = 0; // if there are no wd allocated, it will always be the chunk allocated address

   // Computing Data Device pointers and Data Devicesinfo
   size_DPtrs    = sizeof(DD *) * num_devices;
   offset_DPtrs  = NANOS_ALIGNED_MEMORY_OFFSET(offset_Data, size_Data, __alignof__( DD*) );

   // Computing Copies info
   if ( num_copies != 0 ) {
      size_CopyData = sizeof(CopyData);
      size_Copies   = size_CopyData * num_copies;
      offset_Copies = NANOS_ALIGNED_MEMORY_OFFSET(offset_DPtrs, size_DPtrs, __alignof__(nanos_copy_data_t) );
      // There must be at least 1 dimension entry
      size_Dimensions = num_dimensions * sizeof(nanos_region_dimension_internal_t);
      offset_Dimensions = NANOS_ALIGNED_MEMORY_OFFSET(offset_Copies, size_Copies, __alignof__(nanos_region_dimension_internal_t) );
   } else {
      size_Copies = 0;
      // No dimensions
      size_Dimensions = 0;
      offset_Copies = offset_Dimensions = NANOS_ALIGNED_MEMORY_OFFSET(offset_DPtrs, size_DPtrs, 1);
   }

   // Computing description char * + description
   if ( description == NULL ) {
      offset_DESC = offset_Dimensions;
      size_DESC = size_Dimensions;
   } else {
      offset_DESC = NANOS_ALIGNED_MEMORY_OFFSET(offset_Dimensions, size_Dimensions, __alignof__ (void*) );
      size_DESC = (strlen(description)+1) * sizeof(char);
   }

   // Computing Internal Data info and total size
   static size_t size_PMD   = _pmInterface->getInternalDataSize();
   if ( size_PMD != 0 ) {
      static size_t align_PMD = _pmInterface->getInternalDataAlignment();
      offset_PMD = NANOS_ALIGNED_MEMORY_OFFSET(offset_DESC, size_DESC, align_PMD);
   } else {
      offset_PMD = NANOS_ALIGNED_MEMORY_OFFSET(offset_DESC, size_DESC, 1);
   }

   total_size = NANOS_ALIGNED_MEMORY_OFFSET(offset_PMD, size_PMD, 1);

   chunk = NEW char[total_size];

   // allocating WD and DATA
   if ( *uwd == NULL ) *uwd = (SlicedWD *) chunk;
   if ( outline_data != NULL && *outline_data == NULL ) *outline_data = (chunk + offset_Data);

   // allocating and initializing Device Data pointers
   DD **dev_ptrs = ( DD ** ) (chunk + offset_DPtrs);
   for ( i = 0 ; i < num_devices ; i ++ ) dev_ptrs[i] = ( DD* ) devices[i].factory( devices[i].arg );

   ensure ((num_copies==0 && copies==NULL && num_dimensions==0 && dimensions==NULL) || (num_copies!=0 && copies!=NULL && num_dimensions!=0 && dimensions!=NULL ), "Number of copies and copy data conflict" );

   // allocating copy-ins/copy-outs
   if ( copies != NULL && *copies == NULL ) {
      *copies = ( CopyData * ) (chunk + offset_Copies);
      *dimensions = ( nanos_region_dimension_internal_t * ) ( chunk + offset_Dimensions );
   }

   // Copying description string
   if ( description == NULL ) desc = NULL;
   else desc = (chunk + offset_DESC);

   SlicedWD * wd =  new (*uwd) SlicedWD( *slicer, num_devices, dev_ptrs, outline_data_size, outline_data_align,
                                         outline_data != NULL ? *outline_data : NULL, num_copies, (copies == NULL) ? NULL : *copies, desc );
   // Set WD's socket
   wd->setSocket(  getCurrentSocket() );

   // Set total size
   wd->setTotalSize(total_size );
   
   if ( getCurrentSocket() >= sys.getNumSockets() )
      throw NANOS_INVALID_PARAM;

   // initializing internal data
   if ( size_PMD > 0) {
      _pmInterface->initInternalData( chunk + offset_PMD );
      wd->setInternalData( chunk + offset_PMD );
   }

   // add to workgroup
   if ( uwg != NULL ) {
      WG * wg = ( WG * )uwg;
      wg->addWork( *wd );
   }

   // set properties
   if ( props != NULL ) {
      if ( props->tied ) wd->tied();
      wd->setPriority( dyn_props->priority );
      wd->setFinal ( dyn_props->flags.is_final );
   }
   if ( dyn_props && dyn_props->tie_to ) wd->tieTo( *( BaseThread * )dyn_props->tie_to );
}

/*! \brief Duplicates the whole structure for a given WD
 *
 *  \param [out] uwd is the target addr for the new WD
 *  \param [in] wd is the former WD
 *
 *  \return void
 *
 *  \par Description:
 *
 *  This function duplicates the given WD passed as a parameter copying all the
 *  related data included in the layout (devices ptr, data and DD). First it computes
 *  the size for the layout, then it duplicates each one of the chunks (Data,
 *  Device's pointers, internal data, etc). Finally calls WorkDescriptor constructor
 *  using new and placement.
 *
 *  \sa WorkDescriptor, createWD, createSlicedWD, duplicateSlicedWD
 */
void System::duplicateWD ( WD **uwd, WD *wd)
{
   unsigned int i, num_Devices, num_Copies, num_Dimensions;
   DeviceData **dev_data;
   void *data = NULL;
   char *chunk = 0, *chunk_iter;

   size_t size_CopyData;
   size_t size_Data, offset_Data, size_DPtrs, offset_DPtrs, size_Copies, offset_Copies, size_Dimensions, offset_Dimensions, offset_PMD;
   size_t total_size;

   // WD doesn't need to compute offset, it will always be the chunk allocated address

   // Computing Data info
   size_Data = wd->getDataSize();
   if ( *uwd == NULL ) offset_Data = NANOS_ALIGNED_MEMORY_OFFSET(0, sizeof(WD), wd->getDataAlignment() );
   else offset_Data = 0; // if there are no wd allocated, it will always be the chunk allocated address

   // Computing Data Device pointers and Data Devicesinfo
   num_Devices = wd->getNumDevices();
   dev_data = wd->getDevices();
   size_DPtrs    = sizeof(DD *) * num_Devices;
   offset_DPtrs  = NANOS_ALIGNED_MEMORY_OFFSET(offset_Data, size_Data, __alignof__( DD*) );

   // Computing Copies info
   num_Copies = wd->getNumCopies();
   num_Dimensions = 0;
   for ( i = 0; i < num_Copies; i += 1 ) {
      num_Dimensions += wd->getCopies()[i].getNumDimensions();
   }
   if ( num_Copies != 0 ) {
      size_CopyData = sizeof(CopyData);
      size_Copies   = size_CopyData * num_Copies;
      offset_Copies = NANOS_ALIGNED_MEMORY_OFFSET(offset_DPtrs, size_DPtrs, __alignof__(nanos_copy_data_t) );
      // There must be at least 1 dimension entry
      size_Dimensions = num_Dimensions * sizeof(nanos_region_dimension_internal_t);
      offset_Dimensions = NANOS_ALIGNED_MEMORY_OFFSET(offset_Copies, size_Copies, __alignof__(nanos_region_dimension_internal_t) );
   } else {
      size_Copies = 0;
      // No dimensions
      size_Dimensions = 0;
      offset_Copies = offset_Dimensions = NANOS_ALIGNED_MEMORY_OFFSET(offset_DPtrs, size_DPtrs, 1);
   }

   // Computing Internal Data info and total size
   static size_t size_PMD   = _pmInterface->getInternalDataSize();
   if ( size_PMD != 0 ) {
      static size_t align_PMD = _pmInterface->getInternalDataAlignment();
      offset_PMD = NANOS_ALIGNED_MEMORY_OFFSET(offset_Dimensions, size_Dimensions, align_PMD);
      total_size = NANOS_ALIGNED_MEMORY_OFFSET(offset_PMD,size_PMD,1);
   } else {
      offset_PMD = 0; // needed for a gcc warning
      total_size = NANOS_ALIGNED_MEMORY_OFFSET(offset_Dimensions, size_Dimensions, 1);
   }

   chunk = NEW char[total_size];

   // allocating WD and DATA; if size_Data == 0 data keep the NULL value
   if ( *uwd == NULL ) *uwd = (WD *) chunk;
   if ( size_Data != 0 ) {
      data = chunk + offset_Data;
      memcpy ( data, wd->getData(), size_Data );
   }

   // allocating Device Data
   DD **dev_ptrs = ( DD ** ) (chunk + offset_DPtrs);
   for ( i = 0 ; i < num_Devices; i ++ ) {
      dev_ptrs[i] = dev_data[i]->clone();
   }

   // allocate copy-in/copy-outs
   CopyData *wdCopies = ( CopyData * ) (chunk + offset_Copies);
   chunk_iter = chunk + offset_Copies;
   nanos_region_dimension_internal_t *dimensions = ( nanos_region_dimension_internal_t * ) ( chunk + offset_Dimensions );
   for ( i = 0; i < num_Copies; i++ ) {
      CopyData *wdCopiesCurr = ( CopyData * ) chunk_iter;
      *wdCopiesCurr = wd->getCopies()[i];
      memcpy( dimensions, wd->getCopies()[i].getDimensions(), sizeof( nanos_region_dimension_internal_t ) * wd->getCopies()[i].getNumDimensions() );
      wdCopiesCurr->setDimensions( dimensions );
      dimensions += wd->getCopies()[i].getNumDimensions();
      chunk_iter += size_CopyData;
   }

   // creating new WD 
   //FIXME jbueno (#758) should we have to take into account dimensions?
   new (*uwd) WD( *wd, dev_ptrs, wdCopies, data );

   // Set total size
   (*uwd)->setTotalSize(total_size );
   
   // initializing internal data
   if ( size_PMD != 0) {
      _pmInterface->initInternalData( chunk + offset_PMD );
      (*uwd)->setInternalData( chunk + offset_PMD );
      memcpy ( chunk + offset_PMD, wd->getInternalData(), size_PMD );
   }
}

/*! \brief Duplicates a given SlicedWD
 *
 *  This function duplicates the given as a parameter WD copying all the
 *  related data (devices ptr, data and DD)
 *
 *  \param [out] uwd is the target addr for the new WD
 *  \param [in] wd is the former WD
 */
void System::duplicateSlicedWD ( SlicedWD **uwd, SlicedWD *wd)
{
   unsigned int i, num_Devices, num_Copies;
   DeviceData **dev_data;
   void *data = NULL;
   char *chunk = 0, *chunk_iter;

   size_t size_CopyData;
   size_t size_Data, offset_Data, size_DPtrs, offset_DPtrs;
   size_t size_Copies, offset_Copies, size_PMD, offset_PMD;
   size_t total_size;

   // WD doesn't need to compute offset, it will always be the chunk allocated address

   // Computing Data info
   size_Data = wd->getDataSize();
   if ( *uwd == NULL ) offset_Data = NANOS_ALIGNED_MEMORY_OFFSET(0, sizeof(SlicedWD), wd->getDataAlignment() );
   else offset_Data = 0; // if there are no wd allocated, it will always be the chunk allocated address

   // Computing Data Device pointers and Data Devicesinfo
   num_Devices = wd->getNumDevices();
   dev_data = wd->getDevices();
   size_DPtrs    = sizeof(DD *) * num_Devices;
   offset_DPtrs  = NANOS_ALIGNED_MEMORY_OFFSET(offset_Data, size_Data, __alignof__( DD*) );

   // Computing Copies info
   num_Copies = wd->getNumCopies();
   if ( num_Copies != 0 ) {
      size_CopyData = sizeof(CopyData);
      size_Copies   = size_CopyData * num_Copies;
      offset_Copies = NANOS_ALIGNED_MEMORY_OFFSET(offset_DPtrs, size_DPtrs, __alignof__(nanos_copy_data_t) );
   } else {
      size_Copies = 0;
      offset_Copies = NANOS_ALIGNED_MEMORY_OFFSET(offset_DPtrs, size_DPtrs, 1);
   }

   // Computing Internal Data info and total size
   size_PMD   = _pmInterface->getInternalDataSize();
   if ( size_PMD != 0 ) {
      offset_PMD = NANOS_ALIGNED_MEMORY_OFFSET(offset_Copies, size_Copies, _pmInterface->getInternalDataAlignment());
   } else {
      offset_PMD = NANOS_ALIGNED_MEMORY_OFFSET(offset_Copies, size_Copies, 1);
   }

   total_size = NANOS_ALIGNED_MEMORY_OFFSET(offset_PMD, size_PMD, 1);

   chunk = NEW char[total_size];

   // allocating WD and DATA
   if ( *uwd == NULL ) *uwd = (SlicedWD *) chunk;
   if ( size_Data != 0 ) {
      data = chunk + offset_Data;
      memcpy ( data, wd->getData(), size_Data );
   }

   // allocating Device Data
   DD **dev_ptrs = ( DD ** ) (chunk + offset_DPtrs);
   for ( i = 0 ; i < num_Devices; i ++ ) {
      dev_ptrs[i] = dev_data[i]->clone();
   }

   // allocate copy-in/copy-outs
   CopyData *wdCopies = ( CopyData * ) (chunk + offset_Copies);
   chunk_iter = (chunk + offset_Copies);
   for ( i = 0; i < num_Copies; i++ ) {
      CopyData *wdCopiesCurr = ( CopyData * ) chunk_iter;
      *wdCopiesCurr = wd->getCopies()[i];
      chunk_iter += size_CopyData;
   }

   // creating new SlicedWD 
   new (*uwd) SlicedWD( *(wd->getSlicer()), *((WD *)wd), dev_ptrs, wdCopies, data );

   // Set total size
   (*uwd)->setTotalSize(total_size );
   
   // initializing internal data
   if ( size_PMD != 0) {
      _pmInterface->initInternalData( chunk + offset_PMD );
      (*uwd)->setInternalData( chunk + offset_PMD );
      memcpy ( chunk + offset_PMD, wd->getInternalData(), size_PMD );
   }
}

void System::setupWD ( WD &work, WD *parent )
{
   work.setParent ( parent );
   work.setDepth( parent->getDepth() +1 );
   
   // Inherit priority
   if ( parent != NULL ){
      // Add the specified priority to its parent's
      work.setPriority( work.getPriority() + parent->getPriority() );
   }

   // Prepare private copy structures to use relative addresses
   work.prepareCopies();

   // Invoke pmInterface
   _pmInterface->setupWD(work);
   Scheduler::updateCreateStats(work);
}

void System::submit ( WD &work )
{
   SchedulePolicy* policy = getDefaultSchedulePolicy();
   policy->onSystemSubmit( work, SchedulePolicy::SYS_SUBMIT );
   work.submit();
}

/*! \brief Submit WorkDescriptor to its parent's  dependencies domain
 */
void System::submitWithDependencies (WD& work, size_t numDataAccesses, DataAccess* dataAccesses)
{
   SchedulePolicy* policy = getDefaultSchedulePolicy();
   policy->onSystemSubmit( work, SchedulePolicy::SYS_SUBMIT_WITH_DEPENDENCIES );
   WD *current = myThread->getCurrentWD(); 
   current->submitWithDependencies( work, numDataAccesses , dataAccesses);
}

/*! \brief Wait on the current WorkDescriptor's domain for some dependenices to be satisfied
 */
void System::waitOn( size_t numDataAccesses, DataAccess* dataAccesses )
{
   WD* current = myThread->getCurrentWD();
   current->waitOn( numDataAccesses, dataAccesses );
}


void System::inlineWork ( WD &work )
{
   SchedulePolicy* policy = getDefaultSchedulePolicy();
   policy->onSystemSubmit( work, SchedulePolicy::SYS_INLINE_WORK );
   //! \todo choose actual (active) device...
   if ( Scheduler::checkBasicConstraints( work, *myThread ) ) {
      Scheduler::inlineWork( &work );
   } else {
      Scheduler::submitAndWait( work );
   }
}

void System::createWorker( unsigned p )
{
   NANOS_INSTRUMENT( sys.getInstrumentation()->incrementMaxThreads(); )
   PE *pe = createPE ( "smp", getBindingId( p ), _pes.size() );
   _pes.push_back ( pe );
   BaseThread *thread = &pe->startWorker();
   _workers.push_back( thread );
   ++_targetThreads;

   //Set up internal data
   WD & threadWD = thread->getThreadWD();
   if ( _pmInterface->getInternalDataSize() > 0 ) {
      char *data = NEW char[_pmInterface->getInternalDataSize()];
      _pmInterface->initInternalData( data );
      threadWD.setInternalData( data );
   }
   _pmInterface->setupWD( threadWD );
}

BaseThread * System::getUnassignedWorker ( void )
{
   BaseThread *thread;

   for ( unsigned i = 0; i < _workers.size(); i++ ) {
      thread = _workers[i];
      if ( !thread->hasTeam() || thread->isTaggedToSleep() ) {

         // skip if the thread is not in the mask
         if ( sys.getBinding() && !CPU_ISSET( thread->getCpuId(), &_cpu_active_set) )
            continue;

         // recheck availability with exclusive access
         thread->lock();

         if ( thread->hasTeam() && !thread->isTaggedToSleep() ) {
            // we lost it
            thread->unlock();
            continue;
         }

         thread->reserve(); // set team flag only

         thread->unlock();

         return thread;
      }
   }

   return NULL;
}

BaseThread * System:: getAssignedWorker ( void )
{
   BaseThread *thread;

   ThreadList::reverse_iterator rit;
   for ( rit = _workers.rbegin(); rit != _workers.rend(); ++rit ) {
      thread = *rit;
      if ( thread->hasTeam() && !thread->isTaggedToSleep() ) {

         // recheck availability with exclusive access
         thread->lock();

         if ( !thread->hasTeam() || thread->isTaggedToSleep() ) {
            // we lost it
            thread->unlock();
            continue;
         }

         thread->unlock();

         return thread;
      }
   }

   return NULL;
}

BaseThread * System::getWorker ( unsigned int n )
{
   if ( n < _workers.size() ) return _workers[n];
   else return NULL;
}

void System::acquireWorker ( ThreadTeam * team, BaseThread * thread, bool enter, bool star, bool creator )
{
   int thId = team->addThread( thread, star, creator );
   TeamData *data = NEW TeamData();

   data->setStar(star);

   SchedulePolicy &sched = team->getSchedulePolicy();
   ScheduleThreadData *sthdata = 0;
   if ( sched.getThreadDataSize() > 0 )
      sthdata = sched.createThreadData();

   data->setId(thId);
   data->setTeam(team);
   data->setScheduleData(sthdata);
   if ( creator )
      data->setParentTeamData(thread->getTeamData());

   if ( enter ) thread->enterTeam( data );
   else thread->setNextTeamData( data );

   // The sleep flag must be set before to avoid race conditions
   thread->wakeup();  // set sleep flag only
   thread->reserve(); // set team flag only

   debug( "added thread " << thread << " with id " << toString<int>(thId) << " to " << team );
}

void System::releaseWorker ( BaseThread * thread )
{
   ThreadTeam *team = thread->getTeam();
   unsigned thread_id = thread->getTeamId();

   //! \todo destroy if too many?
   debug("Releasing thread " << thread << " from team " << team );

   if ( _enable_dlb && thread->getTeamId() != 0 ) {
      // teamless threads will only sleep if DLB is enabled
      thread->sleep();
   } else {
      thread->leaveTeam();
      team->removeThread(thread_id);
   }
}

int System::getNumWorkers( DeviceData *arch )
{
   int n = 0;

   for ( ThreadList::iterator it = _workers.begin(); it != _workers.end(); it++ ) {
      if ( arch->isCompatible( ( *it )->runningOn()->getDeviceType() ) ) n++;
   }
   return n;
}

ThreadTeam * System::createTeam ( unsigned nthreads, void *constraints, bool reuseCurrent,
                                  bool enterCurrent, bool enterOthers, bool starringCurrent, bool starringOthers )
{
   if ( nthreads == 0 ) {
      nthreads = getNumThreads();
   }
   
   SchedulePolicy *sched = 0;
   if ( !sched ) sched = sys.getDefaultSchedulePolicy();

   ScheduleTeamData *stdata = 0;
   if ( sched->getTeamDataSize() > 0 )
      stdata = sched->createTeamData();

   // create team
   ThreadTeam * team = NEW ThreadTeam( nthreads, *sched, stdata, *_defBarrFactory(), *(_pmInterface->getThreadTeamData()),
                                       reuseCurrent ? myThread->getTeam() : NULL );

   debug( "Creating team " << team << " of " << nthreads << " threads" );

   // find threads
   if ( reuseCurrent ) {
      acquireWorker( team, myThread, enterCurrent, starringCurrent, /* creator */ true );
      nthreads--;
   }

   while ( nthreads > 0 ) {
      BaseThread *thread = getUnassignedWorker();

      if ( !thread ) {
         createWorker( _pes.size() );
         _numPEs++;
         _numThreads++;
         continue;
      }

      thread->lock();
      acquireWorker( team, thread, enterOthers, starringOthers, /* creator */ false );
      thread->signal();
      thread->unlock();
      nthreads--;
   }

   team->init();

   return team;
}

void System::endTeam ( ThreadTeam *team )
{
   debug("Destroying thread team " << team << " with size " << team->size() );

   dlb_returnCpusIfNeeded();
   while ( team->size ( ) > 0 ) {
      // FIXME: Is it really necessary?
      memoryFence();
   }
   
   fatal_cond( team->size() > 0, "Trying to end a team with running threads");
   
   delete team;
}

void System::updateActiveWorkers ( int nthreads )
{
   NANOS_INSTRUMENT ( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
   NANOS_INSTRUMENT ( static nanos_event_key_t num_threads_key = ID->getEventKey("set-num-threads"); )
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raisePointEvents(1, &num_threads_key, (nanos_event_value_t *) &nthreads); )

   BaseThread *thread;
   ThreadTeam *team = myThread->getTeam();

   /* Increase _numThreads */
   while ( nthreads - _numThreads > 0 ) {

      thread = getUnassignedWorker();
      if ( !thread ) {
         createWorker( _pes.size() );
         _numPEs++;
         continue;
      }

      thread->lock();
      acquireWorker( team, thread, /* enterOthers */ true, /* starringOthers */ false, /* creator */ false );
      thread->signal();
      thread->unlock();
      _numThreads++;
   }

   /* Decrease _numThreads */
   while ( nthreads - _numThreads < 0 ) {
      thread = getAssignedWorker();
      thread->sleep();
      _numThreads--;
   }
}

// Not thread-safe
inline void System::applyCpuMask()
{
   NANOS_INSTRUMENT ( static InstrumentationDictionary *ID = sys.getInstrumentation()->getInstrumentationDictionary(); )
   NANOS_INSTRUMENT ( static nanos_event_key_t num_threads_key = ID->getEventKey("set-num-threads"); )
   NANOS_INSTRUMENT ( nanos_event_value_t num_threads_val = (nanos_event_value_t ) CPU_COUNT(&_cpu_active_set) )
   NANOS_INSTRUMENT ( sys.getInstrumentation()->raisePointEvents(1, &num_threads_key, &num_threads_val); )

   BaseThread *thread;
   ThreadTeam *team = myThread->getTeam();

   for ( unsigned pe_id = 0; pe_id < _pes.size() || _numPEs < (size_t)CPU_COUNT(&_cpu_active_set); pe_id++ ) {

      // Create PE & Worker if it does not exist
      if ( pe_id == _pes.size() ) {
         createWorker( pe_id );
      }

      bool pe_dirty = false;
      int pe_binding = getBindingId( pe_id );
      if ( CPU_ISSET( pe_binding, &_cpu_active_set) ) {

         // This PE should be running
         while ( (thread = _pes[pe_id]->getFirstStoppedThread()) != NULL ) {
            thread->lock();
            acquireWorker( team, thread, /* enterOthers */ true, /* starringOthers */ false, /* creator */ false );
            thread->signal();
            thread->unlock();
            _numThreads++;
            pe_dirty = true;
         }
         if ( pe_dirty ) _numPEs++;

      } else {

         // This PE should not
         while ( (thread = _pes[pe_id]->getFirstRunningThread()) != NULL ) {
            thread->sleep();
            _numThreads--;
            pe_dirty = true;
         }
         if ( pe_dirty ) _numPEs--;
      }
   }
   ensure( _numPEs ==  (size_t)CPU_COUNT(&_cpu_active_set), "applyCpuMask fatal error" );
}

void System::getCpuMask ( cpu_set_t *mask ) const
{
   memcpy( mask, &_cpu_active_set, sizeof(cpu_set_t) );
}

void System::setCpuMask ( const cpu_set_t *mask )
{
   memcpy( &_cpu_active_set, mask, sizeof(cpu_set_t) );
   sys.processCpuMask();
}

void System::addCpuMask ( const cpu_set_t *mask )
{
   CPU_OR( &_cpu_active_set, &_cpu_active_set, mask );
   sys.processCpuMask();
}

inline void System::processCpuMask( void )
{
   // if _bindThreads is enabled, update _bindings adding new elements of _cpu_active_set
   if ( sys.getBinding() ) {
      std::ostringstream oss_cpu_idx;
      oss_cpu_idx << "[";
      for ( int cpu=0; cpu<CPU_SETSIZE; cpu++) {
         if ( CPU_ISSET( cpu, &_cpu_active_set ) ) {

            if ( std::find( _bindings.begin(), _bindings.end(), cpu ) == _bindings.end() ) {
               _bindings.push_back( cpu );
            }

            oss_cpu_idx << cpu << ", ";
         }
      }
      oss_cpu_idx << "]";
      verbose0( "PID[" << getpid() << "]. CPU affinity " << oss_cpu_idx.str() );
      if ( _pmInterface->isMalleable() ) {
         sys.applyCpuMask();
      }
   }
   else {
      verbose0( "PID[" << getpid() << "]. Num threads: " << CPU_COUNT( &_cpu_active_set ) );
      if ( _pmInterface->isMalleable() ) {
         sys.updateActiveWorkers( CPU_COUNT( &_cpu_active_set ) );
      }
   }
}

void System::waitUntilThreadsPaused ()
{
   // Wait until all threads are paused
   _pausedThreadsCond.wait();
}

void System::waitUntilThreadsUnpaused ()
{
   // Wait until all threads are paused
   _unpausedThreadsCond.wait();
}

unsigned System::reservePE ( bool reserveNode, unsigned node, bool & reserved )
{
   // For each available PE
   for ( Bindings::reverse_iterator it = _bindings.rbegin(); it != _bindings.rend(); ++it )
   {
      unsigned pe = *it;
      unsigned currentNode = getNodeOfPE( pe );
      
      // If this PE is in the requested node or we don't need to reserve in
      // a certain node
      if ( currentNode == node || !reserveNode )
      {
         // Ensure there is at least one PE for smp
         if ( _bindings.size() == 1 )
         {
            reserved = false;
            warning( "Cannot reserve PE " << pe << ", there is just one PE left. It will be shared." );
         }
         else
         {
            // Take this pe out of the available bindings list.
            _bindings.erase( --( it.base() ) );
            reserved = true;
         }
         return pe;
      }
   }
   // If we reach this point, there are no PEs available for that node.
   verbose( "reservePE failed for node " << node );
   fatal( "There are no available PEs for the requested node" );
}

void * System::getHwlocTopology ()
{
   return _hwlocTopology;
}

void System::environmentSummary( void )
{
   /* Get Specific Mask String (depending on _bindThreads) */
   cpu_set_t *cpu_set = _bindThreads ? &_cpu_active_set : &_cpu_set;
   std::ostringstream mask;
   mask << "[ ";
   for ( int i=0; i<CPU_SETSIZE; i++ ) {
      if ( CPU_ISSET(i, cpu_set) )
         mask << i << ", ";
   }
   mask << "]";

   /* Get Prog. Model string */
   std::string prog_model;
   switch ( getInitialMode() )
   {
      case POOL:
         prog_model = "OmpSs";
         break;
      case ONE_THREAD:
         prog_model = "OpenMP";
         break;
      default:
         prog_model = "Unknown";
         break;
   }

   message0( "========== Nanos++ Initial Environment Summary ==========" );
   message0( "=== PID:            " << getpid() );
   message0( "=== Num. threads:   " << _numThreads );
   message0( "=== Active CPUs:    " << mask.str() );
   message0( "=== Binding:        " << std::boolalpha << _bindThreads );
   message0( "=== Prog. Model:    " << prog_model );

   for ( ArchitecturePlugins::const_iterator it = _archs.begin();
        it != _archs.end(); ++it ) {

      // Temporarily hide SMP plugin because it has empty information
      if ( strcmp( (*it)->getName(), "SMP PE Plugin" ) == 0 )
         continue;

      message0( "=== Plugin:         " << (*it)->getName() );
      message0( "===  | Threads:     " << (*it)->getNumThreads() );
   }

   message0( "=========================================================" );

   // Get start time
   _summary_start_time = time(NULL);
}

void System::admitCurrentThread ( void )
{
   int pe_id = _pes.size();   

   //! \note Create a new PE and configure it
   PE *pe = createPE ( "smp", getBindingId( pe_id ), pe_id );
   pe->setNUMANode( getNodeOfPE( pe->getId() ) );
   _pes.push_back ( pe );

   //! \note Create a new Thread object and associate it to the current thread
   BaseThread *thread = &pe->associateThisThread ( /* untie */ true ) ;
   _workers.push_back( thread );

   //! \note Update current cpu active set mask
   CPU_SET( getBindingId( pe_id ), &_cpu_active_set );

   //! \note Getting Programming Model interface data
   WD &mainWD = *myThread->getCurrentWD();
   (void) mainWD.getDirectory(true); // FIXME: this may cause future problems
   if ( _pmInterface->getInternalDataSize() > 0 ) {
      char *data = NEW char[_pmInterface->getInternalDataSize()];
      _pmInterface->initInternalData( data );
      mainWD.setInternalData( data );
   }

   //! \note Include thread into main thread
   acquireWorker( _mainTeam, thread, /* enter */ true, /* starring */ false, /* creator */ false );
   
}

void System::expelCurrentThread ( void )
{
   int pe_id =  myThread->runningOn()->getUId();
   _pes.erase( _pes.begin() + pe_id );
   _workers.erase ( _workers.begin() + myThread->getId() );
}

void System::executionSummary( void )
{
//davidp
std::cout << "End of the execution\n";
getStatistics();
//

   time_t seconds = time(NULL) -_summary_start_time;
   message0( "============ Nanos++ Final Execution Summary ============" );
   message0( "=== Application ended in " << seconds << " seconds" );
   message0( "=== " << getCreatedTasks() << " tasks have been executed" );
   message0( "=========================================================" );

	



}

//If someone needs argc and argv, it may be possible, but then a fortran 
//main should be done too
void System::ompss_nanox_main(){
    #ifdef MPI_DEV
    //This function will already do exit(0) after the slave finishes (when we are on slave)
    nanos::ext::MPIProcessor::mpiOffloadSlaveMain();
    #endif    
}
