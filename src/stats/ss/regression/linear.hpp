#ifndef SS_LINEAR_HPP
#define SS_LINEAR_HPP

#include <ss/dist.hpp>
#include <ss/matrix.hpp>
#include <ss/data/errors.hpp>
#include <ss/internal/vargs.hpp>
#include <ss/regression/model.hpp>

namespace SS
{
    struct LinearResults
    {
        std::vector<Coefficient> coeffs;

        Variation error;
        Variation model;
        Variation total;

        // R2
        Real r2;

        // F-statistic for the regression
        Real f;

        // P-value for the regression
        Real p;
    };

    namespace Internal
    {
        /*
         * Computes the least-square fitting by SVD decomposition
         */

        inline LinearResults lSVD(const Matrix &Y, const Matrix &X)
        {
            // The size of the sample
            const auto n = Y.rows();
            
            // Number of coefficients (including the constant coefficient)
            const auto p = X.cols();
            
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
            
            const auto &U = svd.matrixU();
            const auto &V = svd.matrixV();
            
            /*
             * Construct the diagonal matrix for singular values
             */
            
            Matrix D = Eigen::MatrixXd::Constant(X.cols(), X.cols(), 0.0);
            
            for (auto i = 0; i < svd.singularValues().size(); i++)
            {
                D(i,i) = svd.singularValues()[i];
            }
            
            // Calculate for the coefficients
            const Matrix B = V * D.inverse() * U.transpose() * Y;
            
            const Matrix RSE = (Y.transpose() * Y) - (B.transpose() * X.transpose() * Y);
            
            // Degree of freedom for MSE
            const auto errorDF = (n - p);
            
            // Unbiased estimate for the variance
            const auto MSE = RSE(0,0) / errorDF;
            
            // Compute for the variance of the estimators
            const Matrix EV = MSE * (V * (D.inverse() * D.inverse()) * V.transpose());
            
            assert(B.rows() == EV.rows());
            
            LinearResults r;
            
            for (auto i = 0; i < B.rows(); i++)
            {
                Coefficient c;
                
                c.est = B(i, 0);
                c.se  = sqrt(EV(i, i));
                
                /*
                 * Under the null hypothesis, each coefficent has a value of zero. The coefficents
                 * are normally distributed while the population coefficient is unknown.
                 */
                
                c.test = c.se ? c.est / c.se : NAN;
                
                // Each coefficient inherits the same degree of freedom from MSE
                c.df = errorDF;
                
                // Compute p-value for the coefficient
                c.p = !std::isinf(c.test) && !std::isnan(c.test) ? Internal::pval(Internal::pt(c.test, errorDF), TwoSided) : P(NAN);
                
                c.upper = !std::isnan(c.test) ? c.est + c.test * Internal::qt(0.95, errorDF) : NAN;
                c.lower = !std::isnan(c.test) ? c.est - c.test * Internal::qt(0.95, errorDF) : NAN;
                
                r.coeffs.push_back(c);
            }
            
            // Residuals of the the model
            const Matrix E = Y - X * B;
            
            r.error.df = errorDF;
            r.error.ss = (E.transpose() * E)(0, 0);
            r.error.ms = r.error.df ? r.error.ss / r.error.df : NAN;
            
            Real sum = 0.0;
            
            for (auto i = 0; i < Y.col(0).size(); i++)
            {
                sum += (Y.col(0))(i);
            }
            
            const auto mean = sum / Y.col(0).size();
            
            /*
             * Computes for the SSM
             */
            
            #define DEFINE_EY() Matrix EY(Y.rows(), 1); EY.setConstant(mean);
            DEFINE_EY();
            
            const Matrix M = (X * B) - EY;
            
            r.model.df = p - 1;
            r.model.ss = (M.transpose() * M)(0, 0);
            r.model.ms = r.model.ss / r.model.df;
            
            /*
             * Computes for total sum of squares
             */
            
            r.total.df = r.model.df + r.error.df;
            r.total.ss = r.model.ss + r.error.ss;
            r.total.ms = r.total.ss / r.total.df;

            r.f = r.error.ms ? r.model.ms / r.error.ms : NAN;
            r.p = (!std::isinf(r.f) && !std::isnan(r.f) && r.error.df) ? P(1.0 - Internal::pf(r.f, r.model.df, r.error.df)) : P(NAN);

            r.r2 = r.model.ss / r.total.ss;
            
            /*
             * Note that we can't set it earlier because we'd need it for computing r2
             */
            
            if (!r.error.df)
            {
                r.error.ss = NAN;
            }
            
            return r;
        }
    }
    
    inline LinearResults linearModel(const Matrix &Y, const Matrix &X)
    {
        SS_ASSERT(X.rows() >= 2, "Two or more samples is required");
        SS_ASSERT(Y.rows() == X.rows(), "Dimension mismatch in the input variables");

        return Internal::lSVD(Y, X);
    }

    template <typename T, typename... Args> LinearResults linearModel(const T &t, Args... args)
    {
        LinearResults r;
        
        Internal::vArgs([&](const std::vector<const T *> &p)
        {
            /*
             * Create matrix for the dependent variable
             */

            const auto y = *p[0];
            
            Eigen::MatrixXd Y(y.size(), 1);
            Y.col(0) = Eigen::VectorXd::Map(&y[0], y.size());
            
            /*
             * Create matrix for the independent variables
             */
            
            Eigen::MatrixXd X(y.size(), p.size());
            X.setConstant(1.0);

            for (auto i = 1u; i < p.size(); i++)
            {
                const auto x = *p[i];
                X.col(i) = Eigen::VectorXd::Map(&(x[0]), x.size());
            }

            r = linearModel(Y, X);
        }, t, args...);

        return r;
    }
}

#endif
