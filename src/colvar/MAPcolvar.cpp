#include "Colvar.h"
#include "core/Atoms.h"
#include "core/ActionAtomistic.h"
#include "tools/OpenMP.h"
#include "ActionRegister.h"

#include <iostream>
#include <fstream>
#include <omp.h>
#include <string>
#include <deque>
#include "tools/MAPtools.h"

namespace PLMD {
    
    namespace colvar {

        class DIFF: public Colvar {
            
            private:
            
                vector<int> atomsNames;
                
                int nLastFrames, gap;
                
                ForceType forceType;
                
                double forceCoeff, symCoeff, flatCoeff, refSymValue;

                vector<double> refPeaksAngles, refPeaksIntensities, refNullPeaks;
                
                std::deque<vector<Vector>> memForces;
                
                Lattice refLattice;
                Constants constants;
                
                unsigned nt;
                
                struct BoolHandler {
                    bool Push = false;
                    bool Pull = false;
                    bool Sym = false;
                    bool Flat = false;
                    bool isFirstStep = true;
                };
                
                BoolHandler Do;
                
            public:

                static void registerKeywords(Keywords& keys);
                DIFF(const ActionOptions& ao);
                virtual void calculate();
        };

        PLUMED_REGISTER_ACTION(DIFF, "DIFF")

    
        // DEFINE KEYWORDS IN plumed.dat INPUT FILE
        void DIFF::registerKeywords(Keywords& keys) {
            Colvar::registerKeywords(keys);
            
            keys.add("compulsory","NAMES_FILE","Atomic numbers");
            
            keys.add("compulsory","FORCE_COEFF","Force coefficient");
            keys.add("compulsory","FORCE_TYPE","Specifies the type of MD forces calculation (SINGULAR, MEDIAN, MAX)");
            
            keys.addFlag("NO_DIFF",false,"Diffraction algorithm switch");
            keys.add("optional","PATTERN_FILE","Diffraction pattern");
            keys.add("optional","NUM_LAST_FRAMES","Number of the memorized arrays of structures");
            keys.add("optional","GAP","Number of frames between force updates");
            
            keys.add("optional","CELL_A","Reference cell A length (Angstrom)");
            keys.add("optional","CELL_B","Reference cell B length (Angstrom)");
            keys.add("optional","CELL_C","Reference cell C length (Angstrom)");
            keys.add("optional","CELL_ALPHA","Reference cell ALPHA angle (Degrees)");
            keys.add("optional","CELL_BETA","Reference cell BETA angle (Degrees)");
            keys.add("optional","CELL_GAMMA","Reference cell GAMMA angle (Degrees)");
            
            keys.add("optional","LAMBDA_SYM_COEFF","Lasso regression coefficient (Specifies that symmetrization algorithm is active)");    
            
            keys.add("optional","NULL_PEAKS_FILE","Theta angels of the peaks (Specifies that flattening algorithm is active)");                    
            keys.add("optional","LAMBDA_FLAT_COEFF","Flattening coefficient");                    
        }
        
        DIFF::DIFF(const ActionOptions& ao): PLUMED_COLVAR_INIT(ao), nLastFrames(1), gap(1), symCoeff(0.0), flatCoeff(0.0), refSymValue(0.0)
        {
            using std::string;
            using std::to_string;
            using std::stringstream;
            
            std::ifstream in;
            
            int natoms = 0;
            string names, line;
            parse("NAMES_FILE", names);
            in.open(names);
            if(in.is_open()) {
                while(getline(in, line)) {
                    stringstream ss(line);
                    int tempName;
                    ss >> tempName;
                    atomsNames.push_back(tempName - 1);
                    natoms += 1;
                }
                in.close();
            } else { 
                error("Atomic numbers file " + names + " not found"); 
            }
            
            parse("FORCE_COEFF", forceCoeff);
            if(forceCoeff < 0) {
                log.printf("Pulling mode is active \n"); 
                Do.Pull = true;
            } else if(forceCoeff == 0) { 
                error("FORCE_COEFF value cannot be zero"); 
            } else { 
                log.printf("Pushing mode is active \n"); 
                Do.Push = true;
            }
            
            bool tempBool;
            parseFlag("NO_DIFF", tempBool);
            if(tempBool) {
                log.printf("Diffraction mode is unactive \n");
                Do.Push = false;
                Do.Pull = false;
            }
            
            parse("LAMBDA_SYM_COEFF", symCoeff);
            if(symCoeff < 0.0) { 
                error("LAMBDA_SYM_COEFF value cannot be negative"); 
            } else if(symCoeff > 0.0) {
                log.printf("Symmetrization algorithm is active \n");
                Do.Sym = true;
            } 
            
            
            string inForceType;
            parse("FORCE_TYPE", inForceType);
            forceType = ForceType::SINGULAR;
            if(inForceType == "SINGULAR") { forceType = ForceType::SINGULAR; }
            else if (inForceType == "MAX") { forceType = ForceType::MAX; }
            else if (inForceType == "MEDIAN_POINT") { forceType = ForceType::MEDIAN_POINT; }
            else if (inForceType == "MEDIAN_INTEGRAL") { forceType = ForceType::MEDIAN_INTEGRAL; }
            else { error("Unknown force calculation type " + inForceType); }
            
            parse("GAP", gap);
            if(gap <= 0) {
                error("GAP value cannot be negative or zero");
            }
            
            parse("NUM_LAST_FRAMES", nLastFrames);
            if(nLastFrames <= 0) { 
                error("NUM_LAST_FRAMES value cannot be negative or zero");
            }
        

            if(Do.Pull) {
                string pattern;
                parse("PATTERN_FILE", pattern);
                in.open(pattern);
                
                if(in.is_open()) {
                    while(getline(in, line)) {
                        stringstream ss(line);
                        double tempAngle, tempIntensity;
                        ss >> tempAngle >> tempIntensity;
                        refPeaksAngles.push_back(tempAngle);
                        refPeaksIntensities.push_back(tempIntensity);
                    }
                    in.close();
                } else { 
                    error("Pattern file not found");
                }
                
                double a,b,c,alpha,beta,gamma;
                
                parse("CELL_A", a);
                parse("CELL_B", b);
                parse("CELL_C", c);
                parse("CELL_ALPHA", alpha);
                parse("CELL_BETA", beta);
                parse("CELL_GAMMA", gamma);
                
                if(a <= 0 || b <= 0 || c <= 0 || alpha <= 0 || beta <= 0 || gamma <= 0) {
                    error("Reference cell parameters cannot be negative or zero when pulling is active");
                }
                alpha *= DEG_TO_RAD;
                beta *= DEG_TO_RAD;
                gamma *= DEG_TO_RAD;

                refLattice = Lattice::fromParameters(a, b, c, alpha, beta, gamma);
            
                auto maxIntensity = std::max_element(refPeaksIntensities.begin(), refPeaksIntensities.end());
                int npts = static_cast<int>(refPeaksAngles.size());
                
                for (int i = 0; i < npts; ++i) {
                    refPeaksIntensities[i] /= *maxIntensity;
                    refPeaksIntensities[i] *= 100;
                    refSymValue += refPeaksIntensities[i];
                }
                refSymValue *= symCoeff;  
            
                constants.setInitTheta(refPeaksAngles.front() * DEG_TO_RAD);
                constants.setEndTheta(refPeaksAngles.back() * DEG_TO_RAD);
                constants.setNpts(npts);
                
                parse("LAMBDA_FLAT_COEFF", flatCoeff);
                if(flatCoeff < 0.0) { 
                    error("LAMBDA_FLAT_COEFF value cannot be negative"); 
                } else if(flatCoeff > 0.0) {
                    log.printf("Flattening algorithm is active \n");
                    Do.Flat = true;
                }
                
                refNullPeaks = vector<double>(refPeaksAngles.size(), 0.0);
                
                string nullpeaks;
                parse("NULL_PEAKS_FILE", nullpeaks);
                in.open(nullpeaks);
                if(in.is_open()) {
                    while(getline(in, line)) {
                        stringstream ss(line);
                        double temp;
                        int index;
                        ss >> index >> temp;
                        temp /= *maxIntensity;
                        temp *= 100;
                        refNullPeaks[index] = temp;
                    }
                    in.close();
                }
                else if (Do.Flat) { 
                    error("Flattening is active but NULL_PEAKS_FILE file " + nullpeaks + " not found"); 
                }
            }
            
            vector<AtomNumber> atoms;
            for(int i = 1; i < natoms+1; i++){
                AtomNumber d;
                d.setSerial(i);
                atoms.push_back(d);
            }
            requestAtoms(atoms);

            addValueWithDerivatives();
            setNotPeriodic();
            checkRead();
        }

        
        void DIFF::calculate(){
            vector<Vector> atomsCurrent = getPositions();
            Tensor box = getBox();
            box.operator*=(10);
            Tensor intr = inverse(transpose(box));
            for(size_t i = 0; i < atomsCurrent.size(); i++) {
                atomsCurrent[i].operator*=(10);
                atomsCurrent[i] = matmul(intr, atomsCurrent[i]);
            }
            
            Lattice currentLattice;
            
            if(Do.Pull) {
                currentLattice = refLattice;
            } else {
                currentLattice = Lattice::fromVectors(box);
            }
            
            Structure currentStructure(atomsCurrent, atomsNames, currentLattice);
            Powder currentPowder(currentStructure, constants);
            PreProcess pp(currentPowder.getPeaksIntensities(), currentPowder.getPeaksAngles(), constants);
            if(Do.isFirstStep && !Do.Pull) {
                refPeaksIntensities = pp.getPattern();
                Do.isFirstStep = false;
            }
            Correlation correlation(pp.getPattern(), refPeaksIntensities, constants.getStep());  
            
            /*
                Numerical Derivatives
            */
            
            if(Do.isFirstStep && Do.Pull) {
                double hv = 0.000001;
                atomsCurrent[1][0] -= hv;
                Lattice st = Lattice::fromVectors(box);
                Structure cs(atomsCurrent, atomsNames, st);
                Powder pwm(cs, constants);
                PreProcess ppm(pwm.getPeaksIntensities(), pwm.getPeaksAngles(), constants);
                Correlation cm(ppm.getPattern(), refPeaksIntensities, constants.getStep());  
                double dm = cm.getDifference();
                atomsCurrent[1][0] += 2.0*hv;
                Powder pwp(cs, constants);
                PreProcess ppp(pwp.getPeaksIntensities(), pwp.getPeaksAngles(), constants);
                Correlation cp(ppp.getPattern(), refPeaksIntensities, constants.getStep());  
                double dp = cp.getDifference();
                
                vector<Vector> rawDerivatives;
                currentPowder.calcAllDerivativesAtPoint(rawDerivatives, 1);
                pp.calcGridDerivatives(rawDerivatives, true);
                pp.normalizeDerivativesToHundred(rawDerivatives);
                Vector der = correlation.calcCorrelationOverlapDerivatives(rawDerivatives);
                if(std::fabs((dp - dm) / 2.0 / hv - der[0]) > 0.00001) {
                    error("Incorrect derivatives! Check CELL values");
                }
                atomsCurrent[0][0] += hv;
                Do.isFirstStep = false;
            }
        
            setValue(correlation.getDifference());
            
            // double diff = Correlation::pointDifference(pp.getPattern(), refPeaksIntensities);
            // setValue(diff);
            
            vector<Vector> derivatives(atomsCurrent.size(), Vector(0.0, 0.0, 0.0));
            vector<Vector> symDerivatives(atomsCurrent.size(), Vector(0.0, 0.0, 0.0));
            
            #pragma omp parallel for num_threads(OpenMP::getNumThreads())
            for(size_t i = 0; i < atomsCurrent.size(); ++i) {
                vector<Vector> rawDerivatives;
                currentPowder.calcAllDerivativesAtPoint(rawDerivatives, i);
                pp.calcGridDerivatives(rawDerivatives, true);
                pp.normalizeDerivativesToHundred(rawDerivatives);
                if(Do.Pull || Do.Push) {
                    derivatives[i] += correlation.calcCorrelationOverlapDerivatives(rawDerivatives);
                }
                if(Do.Sym) {
                    double symValue = Correlation::symmetrizationValue(pp.getPattern(), symCoeff, 2);
                    Vector symmDer = Correlation::symmetrizationDerivatives(rawDerivatives, pp.getPattern(), symValue, refSymValue, symCoeff, 2);
                    symDerivatives[i] += symmDer;
                }
                // derivatives[i] = Correlation::pointDerivatives(rawDerivatives, pp.getPattern(), refPeaksIntensities);
            }
            
            std::vector<Vector> mdforces;
            PLMD::ActionAtomistic::atoms.getLocalMDForces(mdforces); 
            vector<double> mdf(mdforces.size(), 0.0);
            
            
            switch(forceType) {
                case ForceType::SINGULAR:
                {
                    for(size_t i = 0; i < mdforces.size(); ++i) {
                        mdf[i] = mdforces[i].modulo();
                    }
                    break;
                }
                case ForceType::MAX:
                {
                    auto maxmdf = std::max_element(mdforces.begin(), mdforces.end(), [](Vector va, Vector vb) { return va.modulo() < vb.modulo(); });
                    double maxval = (*maxmdf).modulo();
                    for(size_t i = 0; i < mdforces.size(); ++i) {
                        mdf[i] = maxval;
                    }
                    break;
                }
                case ForceType::MEDIAN_POINT:
                {
                    std::sort(mdforces.begin(), mdforces.end(), [](Vector va, Vector vb) { return va.modulo() < vb.modulo(); });
                    double medianval = mdforces[mdforces.size() / 2].modulo();
                    for(size_t i = 0; i < mdforces.size(); ++i) {
                        mdf[i] = medianval;
                    }
                    break;
                }
                case ForceType::MEDIAN_INTEGRAL:
                {
                    vector<double> mdmodules;
                    double summed = 0.0;
                    for(size_t i = 0; i < mdforces.size(); ++i) {
                        double tempModulo = mdforces[i].modulo();
                        mdmodules.push_back(tempModulo);
                        summed += tempModulo;
                    }
                    std::sort(mdmodules.begin(), mdmodules.end());
                    double tempSum = 0.0;
                    for(size_t i = 0; i < mdmodules.size(); ++i) {
                        tempSum += mdmodules[i];
                        if(tempSum >= summed / 2.0) { 
                            tempSum = mdmodules[i];
                            break;
                        }
                    }
                    for(size_t i = 0; i < mdforces.size(); ++i) {
                        mdf[i] = tempSum;
                    }
                    break;
                }
            }
            
            std::vector<Vector>& f = modifyForces();
            
            if(nLastFrames == 1 && getStep() % gap == 0) {
                for(size_t i = 0; i < atomsCurrent.size(); ++i) {
                    f[i] = mdf[i]*(forceCoeff*derivatives[i] - std::fabs(forceCoeff)*symDerivatives[i]);
                }
                if(Do.Push) {
                    refPeaksIntensities = pp.getPattern();
                    refPeaksAngles = pp.getGrid();
                }
            } else if(getStep() % gap == 0) {
                memForces.push_back(vector<Vector>());
                for(size_t i = 0; i < atomsCurrent.size(); ++i) {
                    memForces.back().push_back(mdf[i]*(forceCoeff*derivatives[i] - std::fabs(forceCoeff)*symDerivatives[i]));
                    Vector result = Vector(0.0, 0.0, 0.0);
                    for(auto arr: memForces) {
                        result += arr[i];
                    }
                    f[i] = result / memForces.size();
                }
                if(static_cast<int>(memForces.size()) == nLastFrames) {
                    memForces.pop_front();
                }
                if(Do.Push) {
                    refPeaksIntensities = pp.getPattern();
                    refPeaksAngles = pp.getGrid();
                }
            }
            
        }
    }
}
