#ifndef NN_ConvolutionalNetwork_hpp
#define NN_ConvolutionalNetwork_hpp

#include <memory>
#include <ostream>

#include <Eigen/Dense>
#include <AffineNetwork.hpp>

namespace NN
{
	class ConvolutionalNetwork
	{
		class Layer
		{
			std::size_t row, col, depth;
			std::unique_ptr<Eigen::Matrix<double, -1, -1>[]> mats;

		public:
			Layer() : row(0), col(0), depth(0), mats(nullptr) {}
			Layer(std::size_t row, std::size_t col, std::size_t depth)
				: row(row), col(col), depth(depth),
				mats(std::make_unique<Eigen::Matrix<double, -1, -1>[]>(depth))
			{
				for (std::size_t i = 0; i < depth; i++)
				{
					mats[i] = Eigen::Matrix<double, -1, -1>(row, col);
				}
			}

			double const &operator()(std::size_t i, std::size_t j, std::size_t k) const { return mats[i](j, k); }
			double &operator()(std::size_t i, std::size_t j, std::size_t k) { return mats[i](j, k); }

			std::size_t const rows() const { return row; }
			std::size_t const cols() const { return col; }
			std::size_t const depth() const { return depth; }

			friend std::ostream &operator<<(std::ostream &os, Layer const &layer)
			{
				if (layer.depth == 0)
				{
					os << "Undefined layer!\n";
					return os;
				}
				os << layer.mats[0];
				for (std::size_t i = 1; i < layer.depth; i++)
				{
					os << "\n\n" << layer.mats[i];
				}
			}
		};
		class Kernel
		{
            std::size_t row, col, in, out;
            std::unique_ptr<double[]> bias;
            std::unique_ptr<std::unique_ptr<Eigen::Matrix<double, -1, -1>[]>[]> mats;

        public:
            Kernel() : row(0), col(0), in(0), out(0), bias(nullptr), mats(nullptr) {}
            Kernel(std::size_t row, std::size_t col, std::size_t in, std::size_t out)
                : row(row), col(col), in(in), out(out),
                bias(std::make_unique<double[]>(out)),
                mats(std::make_unique<std::unique_ptr<Eigen::Matrix<double, -1, -1>[]>[]>(out))
            {
                for (std::size_t i = 0; i < out; i++)
                {
                    mats[i] = std::make_unique<Eigen::Matrix<double, -1, -1>[]>(in);
                    for (std::size_t j = 0; j < in; j++)
                    {
                        mats[i][j] = Eigen::Matrix<double, -1, -1>(row, col);
                    }
                }
            }

            double const &operator()(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const { return mats[i][j](k, l); }
            double &operator()(std::size_t i, std::size_t j, std::size_t k, std::size_t l) { return mats[i][j](k, l); }
            std::size_t const rows() const { return row; }
            std::size_t const cols() const { return col; }
            std::size_t const ins() const { return in; }
            std::size_t const outs() const { return out; }

            friend std::ostream &operator<<(std::ostream &os, Kernel const &kernel)
            {
                for (std::size_t i = 0; i < kernel.outs(); i++)
                {
                    os << "***\n";
                    for (std::size_t j = 0; j < kernel.ins(); j++)
                    {
                        os << kernel.mats[i][j] << '\n';
                    }
                    if (i < kernel.outs() - 1) os << "\n\n";
                }
            }
		};

		std::size_t depth;
		std::unique_ptr<Layer[]> layer;
		std::unique_ptr<Kernel[]> kernel;
		std::unique_ptr<std::size_t[]> zero_padding;

	public:
        /*
            C = number of convolutional layers
            A = number of affine layers
        */
        template<std::size_t C, std::size_t A>
        struct Structure
        {
            std::array<std::array<std::size_t, 2>, C> convolution_size;
            std::array<std::size_t, A> affine_size;
        };
		ConvolutionalNetwork() : depth(0), layer(nullptr), kernel(nullptr), zero_padding(nullptr) {}
		template<std::size_t C, std::size_t A>
		ConvolutionalNetwork(Structure<C, A> const &structure)
			: depth(C + A),
			layer(std::make_unique<Layer[]>(C + 1)),
			kernel(std::make_unique<Kernel[]>(C + 1)),
			zero_padding(std::make_unique<std::size_t[]>(C))
		{
			for (std::size_t i = 0; i < C; i++)
			{
				
			}
		}
	};
}

#endif