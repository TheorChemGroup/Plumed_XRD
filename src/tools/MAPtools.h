#ifndef __PUMED_tools_MAPtools_h
#define __PLUMED_tools_MAPtools_h
#include "MAPconstants.h"
#include "Tensor.h"
#include "Vector.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <numeric>

using std::vector;

namespace PLMD {
    
    enum class ForceType {MAX, MEDIAN_POINT, MEDIAN_INTEGRAL, SINGULAR};
    
    class Constants
    {   
    
        private:
        
            double _lambda = 1.5406;
        
            double _th2ini = 5.0 * DEG_TO_RAD;
            double _th2end = 50.0 * DEG_TO_RAD;
            
            int _npts = 10001;
            
            double _step;
            
            double _sigma = 0.05;
            double _sigma2, _sigmaD;
            
            int _optstep = 200;
            
            void _initStep();
            void _initSigmas();
        
        public:
        
            void setInitTheta(double angle);
            void setEndTheta(double angle);
            void setNpts(int npts);
            void setSigma(double sigma);
            void setLambda(double lambda);
            
            double getInitTheta() const;
            double getEndTheta() const;
            int getNpts() const;
            double getSigma() const;
            double getSigmaSquared() const;
            double getSigmaD() const;
            double getLambda() const;
            double getStep() const;
            int getOptStep() const;
        
            explicit Constants();
        
    };
    
    class Lattice
    {
        private:
        
            double _a, _b, _c, _alpha, _beta, _gamma;

            Tensor _lattice, _reciprocal;
            Vector _diag;      
            
            explicit Lattice(double a, double b, double c
                            , double alpha, double beta, double gamma
                            , const Tensor& lattice)
                        : _a(a), _b(b), _c(c)
                        , _alpha(alpha), _beta(beta), _gamma(gamma)
                        , _lattice(lattice)
                        , _reciprocal(getReciprocalFromLattice(lattice))
                        , _diag(getDiagFromTensor(_reciprocal)) {}
        
        public:
        
            explicit Lattice() {}
            
            const Tensor& getReciprocal() const;
            const Tensor& getLattice() const;
            const Vector& getDiag() const;
            
            static Tensor getReciprocalFromLattice(const Tensor& lattice);
            static Vector getDiagFromTensor(const Tensor& reciprocal);
            
            static Lattice fromVectors(const Tensor& lattice);
            static Lattice fromParameters(double a, double b, double c
                                        , double alpha, double beta, double gamma);
    };
    
    class Structure
    {
        private:
        
            const vector<Vector>& _atoms;
            const vector<int>& _atomsNames;
            const Lattice& _lattice;
           
        public:
        
            explicit Structure(const vector<Vector>& atoms
                             , const vector<int>& atomsNames
                             , const Lattice& lattice)
                         : _atoms(atoms)
                         , _atomsNames(atomsNames)
                         , _lattice(lattice) {}
        
            const vector<Vector>& getAtoms() const;
            const vector<int>& getAtomsNames() const;
            const Lattice& getLattice() const; 
    };
    
    class Powder
    {
        private:
            
            vector<vector<double>> _scatteringFactors;
            vector<vector<double>> _phases;
            vector<double> _peaksIntensities;
            vector<double> _peaksAngles;
            
            vector<Vector> _peaksMillerIndices;
            vector<Vector> _peaksDerivatives;
            
            const Structure& _structure;
            const Constants& _constants = Constants();
            
            int _npeaks;
            
            void generatePeaks();
            
            Vector calcPeakDerivativeAtPoint(int point, int peakNumber);                  
            

        
        public:
            
            explicit Powder(const Structure& structure, const Constants& constants);
        
            void calcAllDerivativesAtPoint(vector<Vector>& emptyVector, int point);
            
            vector<double>& getPeaksIntensities();
            vector<double>& getPeaksAngles();
    };
    
    
    class PreProcess
    {
        private:
        
            const vector<double> _oldPattern;
            const vector<double> _oldGrid;
            const Constants& _constants;
            
            vector<double> _newPattern;
            vector<double> _newNormPattern;
            vector<double> _newGrid;
            
            int _maxPoint;
            double _maxValue;

        public:
        
            explicit PreProcess(const vector<double>& pattern
                              , const vector<double>& grid
                              , const Constants& constants
                              , const vector<double>& newGrid = vector<double>());
        
            void calcGridDerivatives(vector<Vector>& peaksDerivatives, bool optimize);
            
            void normalizeDerivativesToHundred(vector<Vector>& derivatives);
            static void normalizePatternToIntegral(vector<double>& pattern, double step);
            
            const vector<double>& getPattern() const;
            const vector<double>& getUnnormalizedPattern() const;
            const vector<double>& getGrid() const;
            
    };



    class Correlation
    {
        private:
        
            vector<double> _current;
            vector<double> _ref;
            
            double _step;
            
            const vector<double>& _curRaw;
            
            double _cfg, _cgg, _cff, _dfg;
            
            double calcCorrelationFunction(vector<double>& f, vector<double>& g);
            
            Vector calcCorrelationFunctionDerivatives(vector<double>& pattern, const vector<Vector>& tempDers, bool flag);
            
            
        
        public:
            
            explicit Correlation(const vector<double>& current
                               , const vector<double>& ref
                               , double step);
                           
            void calcCorrelationOverlap();
            Vector calcCorrelationOverlapDerivatives(vector<Vector>& gridDerivatives);
            
            double getDifference();
            
            static double pointDifference(const vector<double>& current, const vector<double>& ref);
            static Vector pointDerivatives(vector<Vector>& gridDerivatives, const vector<double>& current, const vector<double>& ref);
            
            static Vector symmetrizationDerivatives(const vector<Vector>& derivatives
												  , const vector<double>& pattern
                                                  , double curValue
                                                  , double refValue
                                                  , double lambda
												  , double jfac);
                                                  
            static double symmetrizationValue(const vector<double>& pattern, double lambda, double jfac); 
    };
    
    
    
}

#endif