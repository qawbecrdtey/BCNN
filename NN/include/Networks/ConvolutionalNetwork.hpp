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
			std::size_t row, col, depth;
			std::unique_ptr<double[]> bias;
			std::unique_ptr<Eigen::Matrix<double, -1, -1>[]> mats;
		public:
			Kernel() : row(0), col(0), depth(0), bias(nullptr), mats(nullptr) {}
			Kernel(std::size_t row, std::size_t col, std::size_t depth)
				: row(row), col(col), depth(depth),
				bias(std::make_unique<double[]>(depth)),
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
			
			friend std::ostream &operator<<(std::ostream &os, Kernel const &kernel)
			{
				if (kernel.depth == 0)
				{
					os << "Undefined kernel!\n";
					return os;
				}
				os << kernel.bias[0] << '\n' << kernel.mats[0];
				for (std::size_t i = 1; i < kernel.depth; i++)
				{
					os << "\n\n" << kernel.bias[i] << '\n' << kernel.mats[i];
				}
				return os;
			}
		};

		std::size_t depth;
		std::unique_ptr<Layer[]> layer;
		std::unique_ptr<Kernel[]> kernel;
		std::unique_ptr<std::size_t[]> zero_padding;

	public:
		ConvolutionalNetwork() : depth(0), layer(nullptr), kernel(nullptr), zero_padding(nullptr) {}
		template<std::size_t C, std::size_t A>
		ConvolutionalNetwork(std::array<std::array<std::size_t, 2>, C> const &convolution_size, std::array<std::size_t, A> const &affine_size)
			: depth(D),
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