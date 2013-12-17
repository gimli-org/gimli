/***************************************************************************
 *   Copyright (C) 2006-2013 by the resistivity.net development team       *
 *   Carsten Rücker carsten@resistivity.net                                *
 *   Thomas Günther thomas@resistivity.net                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "datacontainer.h"
#include "pos.h"
#include "numericbase.h"
#include "vectortemplates.h"

namespace GIMLI{

DataContainer::DataContainer(){
    initDefaults();
    //std::cout << "DataContainer(){" << std::endl;
}

DataContainer::DataContainer(const std::string & fileName, bool sensorIndicesFromOne){
    initDefaults();
    this->load(fileName, sensorIndicesFromOne);
    //std::cout << "DataContainer(const std::string & fileName){" << std::endl;
}

DataContainer::DataContainer(const std::string & fileName, const std::string & sensorTokens, bool sensorIndicesFromOne){
    initDefaults();
    std::vector < std::string > tokenList = getSubstrings(sensorTokens);
    for(Index i=0 ; i < tokenList.size() ; i++) registerSensorIndex(tokenList[i]);
    this->load(fileName, sensorIndicesFromOne);
    //std::cout << "DataContainer(const std::string & fileName){" << std::endl;
}

DataContainer::DataContainer(const DataContainer & data){
    initDefaults();
    this->copy_(data);
    //std::cout << "DataContainer(const DataContainer & data){" << std::endl;
}

DataContainer::~DataContainer(){
    initDefaults();
    clear();
}

DataContainer & DataContainer::operator = (const DataContainer & data){
    if (this != & data){
        this->copy_(data);
    }
    //std::cout << "DataContainer = " << std::endl;
    return * this;
}

void DataContainer::initDefaults(){
    dataMap_["valid"] = RVector(0);
    sensorIndexOnFileFromOne_ = true;
    init();
    initTokenTranslator();
}

void DataContainer::init(){

}

void DataContainer::initTokenTranslator(){
    std::map< std::string, std::string > l;
}

bool DataContainer::haveTranslationForAlias(const std::string & alias) const {
    std::string test(alias);
    test[0] = tolower(test[0]);
    return tT_.count(test) != 0;
}

std::string DataContainer::translateAlias(const std::string & alias) const {
    std::string test(alias);
    test[0] = tolower(test[0]);
    if (tT_.find(test) != tT_.end())
        return tT_.find(test)->second;
    return "no alias";
}

void DataContainer::clear(){
    dataMap_.clear();
}

void DataContainer::copy_(const DataContainer & data){
    clear();

    topoPoints_ = data.additionalPoints();
    sensorPoints_= data.sensorPositions();

    this->resize(data.size());

    inputFormatString_ = data.inputFormatString();
    inputFormatStringSensors_ = data.formatStringSensors();
    dataSensorIdx_  = sensorIdx();

    dataMap_ = data.dataMap();

    dataDescription_ = data.dataDescription();

    tT_ = data.tokenTranslator();

    sensorIndexOnFileFromOne_ = sensorIndexOnFileFromOne();
}

void DataContainer::add(const DataContainer & data, double snap){

    Index start = this->size();
    this->resize(this->size() + data.size());

    IndexArray perm(data.sensorCount(), 0);

    //** merge sensor data
    for (uint i = 0; i < data.sensorCount(); i ++){
        perm[i] = this->createSensor(data.sensorPositions()[i], snap);
    }

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        //** permutate sensorindices
        if (isSensorIndex(it->first)){
            RVector tmp(data.get(it->first));
            for (uint i = 0; i < tmp.size(); i ++){
                SIndex id = (SIndex)tmp[i];
                if (id > -1 && id < (SIndex)perm.size()) {
                    it->second[start + i] = perm[id];
                } else {
                    it->second[start + i] = -1;
                }
            }
        } else {
        //** merge data fields
            it->second.setVal(data.get(it->first), start, -1);
        }
    }
}

long DataContainer::createSensor(const RVector3 & pos, double tolerance){

    long ret = -1;
    for (uint i = 0; i < sensorPoints_.size(); i ++){
        if (pos.distance(sensorPoints_[i]) < tolerance){
            ret = i;
        }
    }

    if (ret == -1){
        ret = sensorPoints_.size();
        sensorPoints_.push_back(pos);
    }
    return ret;
}

void DataContainer::registerSensorIndex(const std::string & token) {
    dataSensorIdx_.insert(token);
    this->set(token, RVector(this->size(), -1.0));
}

bool DataContainer::isSensorIndex(const std::string & token) const {
    return dataSensorIdx_.find(token) != dataSensorIdx_.end();
}

int DataContainer::load(const std::string & fileName, bool sensorIndicesFromOne){
    setSensorIndexOnFileFromOne(sensorIndicesFromOne);
    
	std::fstream file; openInFile(fileName, & file, true);

    std::vector < std::string > row(getNonEmptyRow(file));
	
    if (row.size() != 1){
        throwError(EXIT_DATACONTAINER_NELECS, WHERE_AM_I + " cannot determine data format. " + str(row.size()));
    }

    //** read number of electrodes
    int nSensors = toInt(row[0]);
    if (nSensors < 1){
        throwError(EXIT_DATACONTAINER_NELECS, " cannot determine sensor count " + row[0]);
    }
    RVector x(nSensors, 0.0), y(nSensors, 0.0), z(nSensors, 0.0);

    //** read electrodes format
    //** if no electrodes format is given (no x after comment symbol) take defaults
    std::string sensorFormatDefault("x y z");
    std::vector < std::string > format(getSubstrings(sensorFormatDefault));

    char c; file.get(c);
    if (c == '#') {
        format = getRowSubstrings(file);
//         if (format[0][0] != 'x' && format[0][0] != 'X'){
//             format = getSubstrings(elecsFormatDefault);
//         }
    } else {
        format = getSubstrings(sensorFormatDefault);
    }
    file.unget();

    inputFormatStringSensors_.clear();
    for (uint i = 0; i < format.size(); i ++)
        inputFormatStringSensors_ += format[i] + " ";

    //** read sensor
    for (int i = 0; i < nSensors; i ++){
        row = getNonEmptyRow(file);
        
        if (row.empty()){
            throwError(EXIT_DATACONTAINER_NELECS, 
                       WHERE_AM_I + "To few sensor data. " +
                       str(nSensors) + " Sensors expected but " +
                       str(i) + " found.");
        }

        for (uint j = 0; j < row.size(); j ++){

            if (j == format.size()) break; // no or to few format defined, ignore

            if (format[j] == "x" || format[j] == "X" ||  format[j] == "x/m" || format[j] == "X/m") x[i] = toDouble(row[j]);
            else if (format[j] == "x/mm" || format[j] == "X/mm") x[i] = toDouble(row[j]) / 1000.0;
            else if (format[j] == "y" || format[j] == "Y" ||  format[j] == "y/m" || format[j] == "Y/m" ) y[i] = toDouble(row[j]);
            else if (format[j] == "y/mm" || format[j] == "Y/mm") y[i] = toDouble(row[j]) / 1000.0;
            else if (format[j] == "z" || format[j] == "Z" ||  format[j] == "z/m" || format[j] == "Z/m" ) z[i] = toDouble(row[j]);
            else if (format[j] == "z/mm" || format[j] == "z/mm") z[i] = toDouble(row[j]) / 1000.0;
            else {
                std::cerr << WHERE_AM_I << " Warning! format description unknown: format[" << j << "] = " << format[j] << " column ignored." << std::endl;
                //throwError(EXIT_DATACONTAINER_ELECS_TOKEN, str.str());
            }
        }
    }

    for (int i = 0; i < nSensors; i ++) {
        createSensor(RVector3(x[i], y[i], z[i]).round(1e-12));
    }

    //****************************** Start read the data;
    row = getNonEmptyRow(file);
    if (row.size() != 1) {
        for (Index i = 0; i < row.size(); i ++){
            std::cerr << row[i] << " ";
        }
        std::cerr << std::endl;

        throwError(EXIT_DATACONTAINER_NELECS, WHERE_AM_I + " cannot determine data size. " + str(row.size()));
    }

    int nData = toInt(row[0]);

    //bool fromOne    = true;
    //bool schemeOnly = false;

    if (nData > 0){
        this->resize(nData);

        //** looking for # symbol which start format description section
        file.get(c);
        if (c == '#') {
            format = getRowSubstrings(file);
            if (format.size() == 0){
                throwError(EXIT_DATACONTAINER_NO_DATAFORMAT, WHERE_AM_I + "Can not determine data format.");
            }
//            if (format.size() == 4) schemeOnly = true;
        }
        file.unget();

    }

    std::map< std::string, RVector > tmpMap;

    for (int data = 0; data < nData; data ++){
        row = getNonEmptyRow(file);

        if (row.empty()){
            throwError(EXIT_DATACONTAINER_DATASIZE, 
                       WHERE_AM_I + " To few data. " + str(nData) + 
                       " data expected and " + str(data) + " data found.");
        }

        for (uint j = 0; j < row.size(); j ++){
            if (j == format.size()) break;

            if (!tmpMap.count(format[j])){
                tmpMap.insert(std::pair< std::string, RVector > (format[j], RVector(nData, 0.0))) ;
            }
// 	    std::cout << format[j] << " " << tmpMap[format[j]].size() << std::endl;
// 	    std::cout << row[j] << " " << toDouble(row[j]) << std::endl;
            tmpMap[format[j]][data] = toDouble(row[j]);
        }
    }

    //** renaming formats with the token translator (tT_):
    for (std::map< std::string, RVector >::iterator it = tmpMap.begin();
         it != tmpMap.end(); it++){
        if (haveTranslationForAlias(it->first)){

            double scale = 1.0;
 // should only used in speccialication
            if (it->first.rfind("mV") != std::string::npos ||
                it->first.rfind("mA") != std::string::npos ||
                it->first.rfind("ms") != std::string::npos){
                scale = 1.0 / 1000.0;
            }
            if (it->first.rfind("%") != std::string::npos ){
                scale = 1.0 / 100.0;
            }
            
// 	    std::cout << min(it->second) <<  " " << max(it->second) << std::endl;

//std::cout << "Insert translated:-" << tT_[it->first] << "-" << it->second.size() << std::endl;
            dataMap_[translateAlias(it->first)] = it->second * scale;
        } else {
//std::cout << "Insert:-" <<  it->first << "-" << it->second.size() << std::endl;
            dataMap_[it->first] = it->second;
        }
    }
    
    inputFormatString_.clear();
    for (uint i = 0; i < format.size(); i ++) {
        if (haveTranslationForAlias(format[i])){
            inputFormatString_ += translateAlias(format[i]) + " ";
        } else {
            inputFormatString_ += format[i] + " ";
        }
    }
        
    if (sensorIndexOnFileFromOne_){
        for (std::map< std::string, RVector >::iterator it = dataMap_.begin();
             it!= dataMap_.end(); it ++){
            if (isSensorIndex(it->first)){
                it->second -= RVector(size(), 1.);
            }
        }
    }


    dataMap_["valid"] = 1;

    // validity check should only used in specialization
    this->checkDataValidity();

    //** start read topography;
    row = getNonEmptyRow(file);

    if (row.size() == 1) {
        //** we found topography
        nSensors = toInt(row[0]);
        RVector xt(nSensors), yt(nSensors), zt(nSensors);

        std::string topoFormatDefault("x y z");

        file.get(c);
        if (c == '#') {
            format = getRowSubstrings(file);
            if (format[0] != "x" && format[0] != "X"){
                //** if no electrodes format is given (no x after comment symbol) take defaults;
                format = getSubstrings(topoFormatDefault);
            }
        }
        file.unget();

        //** read topography points;
        for (int i = 0; i < nSensors; i ++){
            row = getNonEmptyRow(file);

            if (row.empty()) {
                throwError(EXIT_DATACONTAINER_NTOPO, WHERE_AM_I 
                           + "To few topo data. " + str(nSensors)
                           + " Topopoints expected and " + str(i) + " found.");
            }

            for (uint j = 0; j < row.size(); j ++){
                if (j == format.size()) break; // no or to few format defined, ignore

                if (     format[j] == "x" || format[j] == "X") xt[i] = toDouble(row[j]);
                else if (format[j] == "y" || format[j] == "Y") yt[i] = toDouble(row[j]);
                else if (format[j] == "z" || format[j] == "Z") zt[i] = toDouble(row[j]);
                else {
                    std::stringstream str;
                    str << " Warning! format description unknown: topo electrode format["
                        << j << "] = " << format[j] << " column ignored." << std::endl;
                    throwError(EXIT_DATACONTAINER_ELECS_TOKEN, str.str());
                }
            }
        }

        for (int i = 0 ; i < nSensors; i ++) {
            topoPoints_.push_back(RVector3(xt[i], yt[i], zt[i]).round(1e-12));
        }
    } // if topo

    file.close();
    return 1;
}

void DataContainer::checkDataValidity(bool remove){
    //** mark inf and nans as invalid
    int nInvalid0 = find(get("valid") < 1).size();
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin();
            it!= dataMap_.end(); it ++){
        this->markInvalid(find(isInfNaN(it->second)));
    }
    int nInvalidNaN = find(get("valid") < 1).size() - nInvalid0;
    if (nInvalidNaN > 0) {
        std::cout << "Warning: removed " << nInvalidNaN << " values due to NaN/Inf!" << std::endl;
    }

    //** check sensor indices < -1 and >= sensorCount()
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
//             std::cout << it->first << " "<< sensorCount()<< " " << find((it->second < -1) | (it->second >= sensorCount())) << std::endl;
//
            this->markInvalid(find((it->second < -1) | (it->second >= sensorCount())));
        }
    }
    int nInvalidIdx = find(get("valid") < 1).size() - nInvalidNaN - nInvalid0;
    if (nInvalidIdx > 0) {
        std::cout << "Warning: removed " << nInvalidIdx << " values due to sensor index out of bounds!" << std::endl;
    }

    //** call local specialization if any
    checkDataValidityLocal();

    if (find(get("valid") < 1).size() > 0){
        std::cout << "Data validity check: found " << find(get("valid") < 1).size() << " invalid data. " << std::endl;
        saveVec(RVector(find(get("valid") < 1)), "invalid.data");
        if (remove) {
            std::cout << "Data validity check: remove invalid data." << std::endl;
            this->removeInvalid();
        }
    }
}

int DataContainer::save(const std::string & fileName, const std::string & formatData,
                        const std::string & formatSensor, bool verbose) const {

    std::fstream file; if (!openOutFile(fileName, & file)) return 0;

    file.precision(14);

    //** START write sensor data
    file << sensorPoints_.size() << std::endl;
    std::string sensorsString(formatSensor);
    file << "# " << sensorsString << std::endl;

    std::vector < std::string > format(getSubstrings(sensorsString));
    for (uint i = 0, imax = sensorPoints_.size(); i < imax; i ++){
        for (uint j = 0; j < format.size(); j ++){
            if (format[j] == "x" || format[j] == "X") file << sensorPoints_[i].x();
            if (format[j] == "y" || format[j] == "Y") file << sensorPoints_[i].y();
            if (format[j] == "z" || format[j] == "Z") file << sensorPoints_[i].z();
            if (j < format.size() -1) file << "\t";
        }
        file << std::endl;
    }

    //** START write data map
    std::vector < const RVector * > outVec;
    std::vector < bool > outInt;

    IndexArray toSaveIdx;
    std::string formatString;

    if (lower(formatData) == "all"){
        toSaveIdx = find(get("valid") > -1);
        formatString = this->tokenList(false);
    } else {
        toSaveIdx = find(get("valid") == 1);
        formatString = formatData;
    }

    int count = toSaveIdx.size();
    file << count << std::endl;
    file << "# " << formatString << std::endl;

    std::vector < std::string > token(getSubstrings(formatString));

    for (uint i = 0; i < token.size(); i ++){
        //  std::cout << token[i] << std::endl;

        std::string valName;
        if (haveTranslationForAlias(token[i])){
            valName = translateAlias(token[i]);
        } else {
            valName = token[i];
        }

        if (dataMap_.count(valName)){
            outVec.push_back(&dataMap_.find(valName)->second );

            if (isSensorIndex(valName)){
                outInt.push_back(true);
            } else {
                outInt.push_back(false);
            }
        } else {
            throwError(1, WHERE_AM_I + " no such data: " + valName);
        }
    }

    for (uint i = 0; i < toSaveIdx.size(); i ++){
        file.setf(std::ios::scientific, std::ios::floatfield);
        file.precision(14);

        for (uint j = 0; j < token.size(); j ++){

            if (outInt[j]){
                file << int((*outVec[j])[toSaveIdx[i]]) + sensorIndexOnFileFromOne_;
            } else {
                 if (token[j] == "valid"){
                    file << int((*outVec[j])[toSaveIdx[i]]);
                 } else {
                    file << (*outVec[j])[toSaveIdx[i]];
                 }
            }
            if (j < token.size() -1) file << "\t";
        }
        file << std::endl;
    }

    //** START write additional points
    file << topoPoints_.size() << std::endl;
    for (uint i = 0; i < topoPoints_.size(); i ++){
        std::cout   << topoPoints_[i].x() << "\t"
                    << topoPoints_[i].y() << "\t"
                    << topoPoints_[i].z() << std::endl;
    }

    file.close();

    if (verbose){
        std::cout << "Wrote: " << fileName << " with " << sensorPoints_.size()
                << " sensors and " << count << " data." <<  std::endl;
    }
    return 1;
}

std::string DataContainer::tokenList(bool withAnnotation) const {
    std::string tokenList;
    if (withAnnotation) tokenList += "SensorIdx: ";
    
    for (std::map< std::string, RVector >::const_iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            tokenList += it->first;
            tokenList += " ";
        }
    }
    if (withAnnotation) tokenList += " Data: ";
    for (std::map< std::string, RVector >::const_iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (!isSensorIndex(it->first)){
            tokenList += it->first;
            tokenList += " ";
        }
    }
    return tokenList;
}

void DataContainer::showInfos() const {
    std::cout << "Sensors: " << this->sensorCount() << ", Data: " << this->size();
    if (topoPoints_.size() > 0){
        std::cout << " Topopoints: " << topoPoints_.size();
    }
    std::cout << std::endl << tokenList() << std::endl;
}

void DataContainer::add(const std::string & token, const RVector & data, const std::string & description){
    this->set(token, data);
    this->setDataDescription(token, description);
//     if (data.size() == this->size()) {
//         dataMap_.insert(make_pair(token, data));
//         this->setDataDescription(token, description);
//     } else {
//         throwError(1, WHERE_AM_I + " wrong data size: " + toStr(this->size()) + " " + toStr(data.size()));
//     }
}

void DataContainer::set(const std::string & token, const RVector & data){
    if (data.size() == this->size()){
        dataMap_[token] = data;
    } else {
        throwError(1, WHERE_AM_I + " wrong data size: " + toStr(this->size()) + " " + toStr(data.size()));
    }
}

const RVector & DataContainer::get(const std::string & token) const {
    if (dataMap_.count(token)) {
        return dataMap_.find(token)->second;
    }

    throwError(1, WHERE_AM_I + " unknown token data for get: " + token + " available are: " + tokenList() );
    return *new RVector(0);
}

const IndexArray DataContainer::id(const std::string & token) const {
    if (!dataMap_.count(token)) {
        throwError(1, WHERE_AM_I + " unknown token data for get: " + token + " available are: " + tokenList() );
    }
        
    if (!isSensorIndex(token)) {
        throwError(1, WHERE_AM_I + " token: " + token + " is not an index list: " + tokenList() );
    }
        
    //RVector idx = dataMap_.find(token)->second;
    IndexArray ret(dataMap_.find(token)->second.size());
    
    for (uint i = 0; i < ret.size(); i ++){
        ret[i] = Index(dataMap_.find(token)->second[i]);
    }
            
    return ret;
}

RVector * DataContainer::ref(const std::string & token){
    if (dataMap_.count(token)) {
        return &dataMap_.find(token)->second;
    }

    throwError(1, WHERE_AM_I + " unknown token data for ref: " + token + " available are: " + tokenList() );
    return NULL;
}

void DataContainer::setDataDescription(const std::string & token, const std::string & description){
//    std::cout << "(this->exists(token)) " << this->exists(token) << std::endl;
    if (this->exists(token)){
        dataDescription_[token] = description;
    }
}

std::string DataContainer::dataDescription(const std::string & token) const {
    if (this->exists(token) && (dataDescription_.count(token))) {
        return dataDescription_.find(token)->second;
    }
    return "";
}

void DataContainer::resize(uint size) {
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin();
            it!= dataMap_.end(); it ++){

        if (isSensorIndex(it->first)){
            // if the data field represent a sensorIDX fill with -1
            it->second.resize(size, -1.0);
        } else {
            // else pure data, fill with 0.0
            it->second.resize(size, 0.0);
        }

    }
}

void DataContainer::removeInvalid(){
    std::vector< size_t > validIdx(find(get("valid") == 1));

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        dataMap_[it->first] = it->second(validIdx);
    }
}

void DataContainer::remove(const IndexArray & idx){
    this->markValid(idx, false);
    this->removeInvalid();
}

DataContainer DataContainer::filter(const IndexArray & idx) const {
    DataContainer data(*this);
    data.markValid(find(data("valid") > -1), false);
    data.markValid(idx, true);
    data.removeInvalid();
    return data;
}

// START Sensor related section
void DataContainer::removeSensorIdx(uint idx){
    IndexArray i(1, idx);
    this->removeSensorIdx(i);
}

void DataContainer::removeSensorIdx(const IndexArray & idx){
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            for (IndexArray::const_iterator id = idx.begin(); id != idx.end(); id ++){
                this->markValid(find(it->second == *id), false);
            }
        }
    }
    this->removeInvalid();
    this->removeUnusedSensors();
}

void DataContainer::removeUnusedSensors(bool verbose){

    uint oldSensorCount = this->sensorCount();

    BVector activeSensors(this->sensorCount(), false);

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            for (uint i = 0; i < it->second.size(); i ++){
                SIndex id = (SIndex)it->second[i];
                if (id > -1 && id < (SIndex)this->sensorCount()) activeSensors[id] = true;
            }
        }
    }

    IndexArray perm(this->sensorCount(), 0);

    std::vector < RVector3 > tmp(sensorPoints_);
    sensorPoints_.clear();

    for (size_t i = 0; i < activeSensors.size(); i ++) {
        if (activeSensors[i]){
            sensorPoints_.push_back(tmp[i]);
        }
        perm[i] = sensorPoints_.size() -1;
    }

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            for (uint i = 0; i < it->second.size(); i ++){
                SIndex id = (SIndex)it->second[i];
                if (id > -1 && id < (SIndex)perm.size()) it->second[i] = perm[id];
            }
        }
    }

    if (verbose){
        std::cout << "Removed " << oldSensorCount - this->sensorCount() << " sensors." << std::endl;
    }
}

bool idPosLesserX(const std::pair < RVector3, Index > & a, const std::pair < RVector3, Index > & b){
    return posLesserX(a.first, b.first);
}

void DataContainer::sortSensorsX(){
    std::vector < std::pair < RVector3, Index > > permSens(sensorCount());
    for (uint i = 0; i < sensorCount(); i ++) permSens[i] = std::pair< RVector3, Index >(sensorPoints_[i], i);

    std::sort(permSens.begin(), permSens.end(), idPosLesserX);

    IndexArray perm(sensorCount());
    for (uint i = 0; i < perm.size(); i ++){
        sensorPoints_[i] = permSens[i].first;
        perm[permSens[i].second] = i;
    }

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            for (uint i = 0; i < it->second.size(); i ++){
                SIndex id = (SIndex)it->second[i];
                if (id > -1 && id < (SIndex)perm.size()) it->second[i] = perm[id];
            }
        }
    }
}

bool ididLesser(const std::pair < Index, Index > & a, const std::pair < Index, Index > & b){
    return a.first < b.first;
}

void DataContainer::sortSensorsIndex(){

    std::vector < std::pair < Index, Index > > permSens(this->size());
    Index nSensorsIdx = dataSensorIdx_.size();

    for (uint i = 0; i < this->size(); i ++) {
        Index sensorUniqueID = 0;
        Index count = nSensorsIdx;
        for (std::set< std::string >::iterator it = dataSensorIdx_.begin(); it!= dataSensorIdx_.end(); it ++){
            count --;
            sensorUniqueID += ((Index)dataMap_[*it][i] + 1) * (Index)powInt(this->sensorCount(), count);
        }
        permSens[i] = std::pair< Index, Index >(sensorUniqueID, i);
    }

    std::sort(permSens.begin(), permSens.end(), ididLesser);

    IndexArray perm(this->size());
    for (uint i = 0; i < perm.size(); i ++){
        perm[i] = permSens[i].second ;
    }

    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        it->second = it->second(perm);
    }
}

void DataContainer::markInvalidSensorIndices(){
    for (std::map< std::string, RVector >::iterator it = dataMap_.begin(); it!= dataMap_.end(); it ++){
        if (isSensorIndex(it->first)){
            this->markValid(find(it->second >= this->sensorCount()), false);
        }
    }
}

void DataContainer::translate(const RVector3 & trans){
    for (uint i = 0; i < sensorPoints_.size(); i ++) {
        sensorPoints_[i].translate(trans);
    }
}

void DataContainer::scale(const RVector3 & scale){
    for (uint i = 0; i < sensorPoints_.size(); i ++) {
        sensorPoints_[i].scale(scale);
    }
}

// END Sensor related section

} // namespace GIMLI{

