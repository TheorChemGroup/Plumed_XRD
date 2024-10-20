#include "MAPtools.h"
#include <iostream>

using std::vector;
using std::sqrt;
using std::sin;
using std::cos;

namespace PLMD {
    
    /*
        CLASS CONSTANTS
    */
    
    Constants::Constants()
    {
        _initStep();
        _initSigmas();
    }
    
    void Constants::_initStep() 
    {
        _step = (_th2end-_th2ini) / (_npts-1) / DEG_TO_RAD;
    }
    
    void Constants::_initSigmas() 
    {
        _sigma2 = _sigma * _sigma;
        _sigmaD = 1.0 / (_sigma2 * 2.0);
    }
    
    void Constants::setInitTheta(double angle) 
    {
        _th2ini = angle;
        _initStep();
    }
    
    void Constants::setEndTheta(double angle)
    {
        _th2end = angle;
        _initStep();
    }
    
    void Constants::setNpts(int npts)
    {
        _npts = npts;
        _initStep();
    }
    
    void Constants::setSigma(double sigma)
    {
        _sigma = sigma;
        _initSigmas();
    }
    
    void Constants::setLambda(double lambda)
    {
        _lambda = lambda;
    }
    
    double Constants::getInitTheta() const { return _th2ini; }
    double Constants::getEndTheta() const { return _th2end; }
    int Constants::getNpts() const { return _npts; }
    double Constants::getSigma() const { return _sigma; }
    double Constants::getSigmaSquared() const { return _sigma2; }
    double Constants::getSigmaD() const { return _sigmaD; }
    double Constants::getLambda() const { return _lambda; }
    double Constants::getStep() const { return _step; }
    int Constants::getOptStep() const { return _optstep; }

    
    /*
        CLASS LATTICE
    */
    
    Lattice Lattice::fromVectors(const Tensor& lattice)
    {
        double asq = 0.0, bsq = 0.0, csq = 0.0;
        double dotAB = 0.0, dotAC = 0.0, dotBC = 0.0;
        for(int i = 0; i < 3; ++i) {
            asq += lattice[0][i] * lattice[0][i];
            bsq += lattice[1][i] * lattice[1][i];
            csq += lattice[2][i] * lattice[2][i];
            dotAB += lattice[0][i] * lattice[1][i];
            dotAC += lattice[0][i] * lattice[2][i];
            dotBC += lattice[1][i] * lattice[2][i];
        }
        
        double a = sqrt(asq);
        double b = sqrt(bsq);
        double c = sqrt(csq);

        using std::acos;
        double alpha = acos(dotBC/bsq/csq);
        double beta = acos(dotAC/asq/csq);
        double gamma = acos(dotAB/asq/bsq);
                       
        return Lattice(a, b, c, alpha, beta, gamma, lattice);
    }

    Lattice Lattice::fromParameters(double a, double b, double c, double alpha, double beta, double gamma)
    {
        double cx = c*cos(beta);
        double cy = c*((cos(alpha) - cos(beta)*cos(gamma))/sin(gamma));
        Tensor tempTensor = Tensor{a, 0.0, 0.0, b*cos(gamma), b*sin(gamma), 0.0, cx, cy, sqrt(c*c - cx*cx - cy*cy)};
        
        return Lattice(a, b, c, alpha, beta, gamma, tempTensor);
    }
    
    Tensor Lattice::getReciprocalFromLattice(const Tensor& lattice) 
    {
        Tensor reciprocal = matmul(lattice, transpose(lattice));
        reciprocal = reciprocal.inverse();
        reciprocal.operator*=(BOHR_TO_A*BOHR_TO_A);
        return reciprocal;
    }
    
    Vector Lattice::getDiagFromTensor(const Tensor& reciprocal)
    {
        Vector diag;
        for(int i = 0; i < 3; ++i) {
            diag[i] = sqrt(reciprocal[i][i]);
        }
        return diag;
    }

    const Vector& Lattice::getDiag() const { return _diag; }
    const Tensor& Lattice::getReciprocal() const { return _reciprocal; }
    const Tensor& Lattice::getLattice() const { return _lattice; }


    /*
        CLASS STRUCTURE
    */
    
    const vector<Vector>& Structure::getAtoms() const { return _atoms; }
    const vector<int>& Structure::getAtomsNames() const { return _atomsNames; }
    const Lattice& Structure::getLattice() const { return _lattice; }

    
    /*
        CLASS POWDER
    */
    
    Powder::Powder(const Structure& structure
                 , const Constants& constants): _structure(structure)
                                              , _constants(constants)
    {
        generatePeaks();
    }

    void Powder::generatePeaks() 
    {
        using std::fabs;
        
        double ieps = 1e-5;
        double iepscont = 1e-10;
        double intmax = 1e15;
        
        double tshift = _constants.getSigma() * sqrt(fabs(-2.0 * std::log(iepscont / intmax))) * DEG_TO_RAD;
        double lambdaBohr = _constants.getLambda() / BOHR_TO_A;

        double smax = sin((_constants.getEndTheta() + tshift) / 2.0);
        
        const Vector& diag = _structure.getLattice().getDiag();
        int hmax = 2 * (int) ceil(2.0 * smax / lambdaBohr / std::min({diag[0], diag[1], diag[2]}));

        _npeaks = 0;
        for (int hcell = 1; hcell <= hmax; ++hcell) {
            for (int h = -hcell; h <= hcell; ++h) {
                for (int k = -hcell; k <= hcell; ++k) {
                    for (int l = -hcell; l <= hcell; ++l) {
                        if (abs(h) != hcell && abs(k) != hcell && abs(l) != hcell) {continue;}

                        Vector hvec = {(double) h, (double) k, (double) l};
                        double dh2 = dotProduct(matmul(_structure.getLattice().getReciprocal(), hvec), hvec);
                        double dh = sqrt(dh2);
                        double dh3 = dh2 * dh;
                        
                        double sth = 0.5 * lambdaBohr * dh;
                        if (fabs(sth) > smax) {continue;}

                        double th = std::asin(sth);
                        double th2 = 2.0 * th;

                        if (th2 < (_constants.getInitTheta() - tshift) || th2 > (_constants.getEndTheta() + tshift)) {continue;}

                        double sthlam = dh / BOHR_TO_A / 2.0;

                        Vector kvec = hvec * M_2PI;

                        double cterm = 0.0;
                        double sterm = 0.0;
                        
                        const vector<int>& atomsNames = _structure.getAtomsNames();
                        const vector<Vector>& atoms = _structure.getAtoms();
                        vector<double> temp_phases, temp_factors;
                        
                        for (size_t i = 0; i < atoms.size(); ++i) {
                            double ffac = 0.0;
                            if (dh < 2) {
                                for(int j = 0; j < 8; j += 2) {
                                    ffac += SFACTOR_1[atomsNames[i]][j] * exp(-SFACTOR_1[atomsNames[i]][j+1] * dh2);
                                }
                                ffac += SFACTOR_1[atomsNames[i]][8];
                            } else {
                                ffac = exp(SFACTOR_2[atomsNames[i]][0] + 
                                           SFACTOR_2[atomsNames[i]][1] * dh +
                                           SFACTOR_2[atomsNames[i]][2] * dh2 / 10.0 +
                                           SFACTOR_2[atomsNames[i]][3] * dh3 / 100.0);
                            }

                            ffac *= exp(-(sthlam * sthlam));
                            
                            temp_factors.push_back(ffac);

                            double phase = dotProduct(kvec, atoms[i]);
                            temp_phases.push_back(phase);
                            
                            cterm += ffac * cos(phase);
                            sterm += ffac * sin(phase);
                        }

                        double res = (cterm * cterm + sterm * sterm) * (0.75 + 0.25 * cos(2.0 * th2)) / sin(th2) / sin(th);

                        if (res > ieps) {
                            _npeaks += 1;
                            _peaksIntensities.push_back(res);
                            _peaksAngles.push_back(th2);
                            _peaksMillerIndices.push_back(hvec);
                            _scatteringFactors.push_back(temp_factors);
                            _phases.push_back(temp_phases);
                        }
                    }
                }
            }
        }
    }
            
    Vector Powder::calcPeakDerivativeAtPoint(int point, int peakNumber)
    {
        double realpart = 0.0;
        double impart = 0.0;
        
        for (size_t i = 0; i < _structure.getAtoms().size(); ++i) {
            realpart += _scatteringFactors[peakNumber][i] * cos(_phases[peakNumber][i]);
            impart += _scatteringFactors[peakNumber][i] * sin(_phases[peakNumber][i]);
        }
        
        double lpf = (1.5 + 0.5 * cos(2.0 * _peaksAngles[peakNumber])) / (sin(_peaksAngles[peakNumber]) * sin(_peaksAngles[peakNumber] / 2.0));
        double tempPart = M_2PI * _scatteringFactors[peakNumber][point] * (-sin(_phases[peakNumber][point]) * realpart +
                                                                               cos(_phases[peakNumber][point]) * impart);
        
        return _peaksMillerIndices[peakNumber] * tempPart * lpf;
    }
    
    void Powder::calcAllDerivativesAtPoint(vector<Vector>& emptyVector, int point) 
    {
        emptyVector.reserve(_peaksIntensities.size());
        for(size_t i = 0; i < _peaksIntensities.size(); ++i) {
            emptyVector.push_back(calcPeakDerivativeAtPoint(point, i));
        }
    }
    
    vector<double>& Powder::getPeaksIntensities() 
    {
        return _peaksIntensities;
    }
    vector<double>& Powder::getPeaksAngles()
    {
        return _peaksAngles;
    }
    
    
    /*
        CLASS PREPROCESS
    */
    
    PreProcess::PreProcess(const vector<double>& pattern
                         , const vector<double>& grid
                         , const Constants& constants
                         , const vector<double>& newGrid)
                     : _oldPattern(pattern)
                     , _oldGrid(grid)
                     , _constants(constants) 
    {
        if(!newGrid.empty()) {
            _newGrid = newGrid;
        }
        else {
            for(int i = 0; i < _constants.getNpts(); ++i) {
                _newGrid.push_back(_constants.getStep()*i + _constants.getInitTheta()/DEG_TO_RAD);
            }
        }
        _newPattern = vector<double>(_newGrid.size(), 0.0);
        for(size_t j = 0; j < _oldGrid.size(); ++j) {
            double tempAngle = _oldGrid[j] / DEG_TO_RAD;
            for(size_t i = 0; i < _newGrid.size(); ++i) {
                _newPattern[i] += _oldPattern[j] * std::exp(-std::pow(_newGrid[i] - tempAngle, 2.0) * _constants.getSigmaD());
            }
        }
        double tempMax = -1.0;
        for(size_t i = 0; i < _newPattern.size(); ++i) {
            if(_newPattern[i] > tempMax) {
                tempMax = _newPattern[i];
                _maxPoint = i;
            }
        }
        _newNormPattern = _newPattern;
        for(size_t i = 0; i < _newPattern.size(); ++i) {
            _newNormPattern[i] *= 100.0 / tempMax;
        }
    }
    
    void PreProcess::calcGridDerivatives(vector<Vector>& derivatives, bool optimize) 
    {
        vector<Vector> newDerivatives(_newGrid.size(), {0.0, 0.0, 0.0});
        if(optimize) {
            for (size_t k = 0; k < derivatives.size(); ++k) {
                double angle = _oldGrid[k] / DEG_TO_RAD;
                int start = static_cast<int>(std::fmax(0.0, floor((angle - _constants.getInitTheta()/DEG_TO_RAD)/_constants.getStep()) - _constants.getOptStep()/2));
                int end = static_cast<int>(std::fmin(start + _constants.getOptStep(), _constants.getNpts()));
                for (int i = start; i < end; ++i) {
                    double co = std::exp(_constants.getSigmaD() * -std::pow(angle-_newGrid[i], 2.0));
                    newDerivatives[i] += derivatives[k]*co;
                }
            }
        }
        else {
            for(size_t k = 0; k < derivatives.size(); ++k) {
                for(size_t i = 0; i < _newGrid.size(); ++i) {
                    newDerivatives[i] += derivatives[k] * std::exp(_constants.getSigmaD() * -std::pow(_newGrid[i] - _oldGrid[k] / DEG_TO_RAD, 2.0));
                }
            }
        }
        derivatives = std::move(newDerivatives);
    }
    
    void PreProcess::normalizeDerivativesToHundred(vector<Vector>& derivatives) 
    {
        Vector tempMax = derivatives[_maxPoint];
        double maxSq = _newPattern[_maxPoint]*_newPattern[_maxPoint];
        for(size_t i = 0; i < derivatives.size(); ++i) {
            derivatives[i] = 100.0 * (derivatives[i]*_newPattern[_maxPoint] - _newPattern[i]*tempMax)/maxSq;
        }
    }
    
    void PreProcess::normalizePatternToIntegral(vector<double>& pattern, double step) 
    {
        double tini = pattern.front() * pattern.front();
        double tend = pattern.back() * pattern.back();
        double summ0 = 0.0;
        for (size_t i = 2; i <= pattern.size()-1; ++i) {
            summ0 += pattern[i - 1] * pattern[i - 1];
        }
        double nor = sqrt((2.0 * summ0 + tini + tend) * step/2.0);
        for (size_t i = 0; i < pattern.size(); ++i) {
            pattern[i] /= nor;
        }
    }
    
    const vector<double>& PreProcess::getPattern() const 
    {
        return _newNormPattern;
    }
    
    const vector<double>& PreProcess::getUnnormalizedPattern() const 
    {
        return _newPattern;
    }
    
    const vector<double>& PreProcess::getGrid() const 
    {
        return _newGrid;
    }
   
    
    /*
        CLASS CORRELATION
    */
    
    Correlation::Correlation(const vector<double>& current
                   , const vector<double>& ref
                   , double step)
               : _current(current)
               , _ref(ref)
               , _step(step) 
               , _curRaw(current)
    {
        PreProcess::normalizePatternToIntegral(_current, _step);
        PreProcess::normalizePatternToIntegral(_ref, _step);
        calcCorrelationOverlap();
    }
        

    double Correlation::calcCorrelationFunction(vector<double>& f, vector<double>& g) 
    {
        double summ = 0.0;
        double w = 0.0;

        for (int i = 1; i <= static_cast<int>(std::floor(1.0 / _step)); ++i) {
            w = std::fmax(1 - static_cast<double>(i * _step), 0.0);
            for (size_t k = 1; k <= f.size() - i; ++k) {
                summ += w * (f[k - 1] * g[k + i - 1] + g[k - 1] * f[k + i - 1]);
            }
        }
        for (size_t i = 1; i <= f.size(); ++i) {
            summ += f[i - 1] * g[i - 1];
        }

        return summ * _step * _step;
    }
    

    void Correlation::calcCorrelationOverlap() 
    {
        _cfg = calcCorrelationFunction(_current, _ref);
        _cgg = calcCorrelationFunction(_current, _current);
        _cff = calcCorrelationFunction(_ref, _ref);
        _dfg = std::fmax(1.0 - _cfg / sqrt(_cff * _cgg), 0.0);
    }
    
    Vector Correlation::calcCorrelationFunctionDerivatives(vector<double>& pattern, const vector<Vector>& tempDers, bool flag) 
    {
        Vector ders;
        for (int i = 1; i <= static_cast<int>(std::floor(1.0 / _step)); ++i) {
            double omega = std::max(1.0 - i * _step, 0.0);
            for (size_t k = 1; k <= pattern.size() - i; ++k) {
                for(int j = 0; j < 3; ++j) {
                    ders[j] += omega * (pattern[k - 1] * tempDers[k + i - 1][j] + tempDers[k - 1][j] * pattern[k + i - 1]);
                }
            }
        }
        for (size_t i = 1; i <= pattern.size(); ++i) {
            for(int j = 0; j < 3; ++j) {
                ders[j] += (pattern[i - 1] * tempDers[i - 1][j]);
            }
        }
        double hhflag = _step * _step * (flag + 1);
        for(int i = 0; i < 3; ++i) {
            ders[i] *= hhflag;
        }

        return ders;
    }
    
    double Correlation::getDifference() 
    {
        return _dfg;
    }
    
    
    Vector Correlation::calcCorrelationOverlapDerivatives(vector<Vector>& preprocessedDers) 
    {
        double tempSum = 0.0;
        for (size_t i = 1; i < _curRaw.size()-1; ++i) {
            tempSum += _curRaw[i] * _curRaw[i];
        }
        double nor = (2.0 * tempSum + _curRaw.front() * _curRaw.front() + _curRaw.back() * _curRaw.back()) * _step/2.0;        
        double norDer[3] = {0.0, 0.0, 0.0};
        for (size_t i = 1; i < _curRaw.size()-1; ++i) {
            for(int j = 0; j < 3; ++j) {
                norDer[j] += preprocessedDers[i][j] * _curRaw[i];
            }
        }
        for(int j = 0; j < 3; ++j) {
            norDer[j] *= 2;
            norDer[j] += preprocessedDers.front()[j] * _curRaw.front() + preprocessedDers.back()[j] * _curRaw.back();
            norDer[j] *= _step;
        }
        vector<Vector> tempDers;

        double sq = sqrt(nor);

        for (size_t i = 0; i < _curRaw.size(); ++i) {
            Vector tempVec;
            for(int j = 0; j < 3; ++j) {
                tempVec[j] = (preprocessedDers[i][j] * sq - (norDer[j] * _curRaw[i]) / (2.0 * sq)) / nor;
            }
            tempDers.push_back(tempVec);
        }
        Vector cfgDers = calcCorrelationFunctionDerivatives(_ref, tempDers, false);
        Vector cggDers = calcCorrelationFunctionDerivatives(_current, tempDers, true);
        
        double sj = sqrt(_cff * _cgg);
        double sg = sqrt(_cff / _cgg);
        
        Vector ders;
        for(int i = 0; i < 3; ++i) {
            ders[i] = (0.5 * _cfg * cggDers[i] * sg - cfgDers[i] * sj) / (_cff * _cgg);
        }
        return ders;
    }
    
    double Correlation::pointDifference(const vector<double>& current, const vector<double>& ref)
    {
        double result = 0.0;
        for(size_t i = 0; i < current.size(); ++i) {
            result += std::pow(current[i] - ref[i], 2.0);
        }
        return result;
        
    }
    
    Vector Correlation::pointDerivatives(vector<Vector>& gridDerivatives, const vector<double>& current, const vector<double>& ref)
    {
        Vector result = Vector(0.0, 0.0, 0.0);
        for(size_t i = 0; i < current.size(); ++i) {
            result += 2.0 * (current[i] - ref[i]) * gridDerivatives[i];
        }
        return result;
    }
    
    Vector Correlation::symmetrizationDerivatives(const vector<Vector>& derivatives, const vector<double>& pattern, double curValue, double refValue, double lambda, double jfac)
    {   
        Vector result = Vector(0.0, 0.0, 0.0);
		double jfac_sq = jfac * jfac;
		
		for(size_t i = 0; i < pattern.size(); ++i) {
			result += derivatives[i] * jfac / (1.0 + jfac_sq * pattern[i] * pattern[i]);
		}
		
        result *= lambda;
        
        if(curValue - refValue < 0) { result *= -1; }
        return result;
    }
    
    double Correlation::symmetrizationValue(const vector<double>& pattern, double lambda, double jfac) 
    {
		double value = 0.0;
		
        for(double intensity: pattern) {
			value += std::atan(jfac * intensity);
		}
        
        value *= lambda;
        
        return value;
    }
    
}
