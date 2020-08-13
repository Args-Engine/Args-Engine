#pragma once
#include <core/types/primitives.hpp>


namespace args::core::math
{
	enum class precision
	{
		bit32,
		bit64,
		lots,
		whole_only
	};


	template <precision p>
	struct precision_chooser
	{
		using type = 
		/*if*/		std::conditional_t<p == precision::bit32,f32,
		/*elseif*/	std::conditional_t<p == precision::bit64,f64,
		/*elseif*/	std::conditional_t<p == precision::lots, f80,
					int32>>>;
	};

	template <precision p>
	using precision_chooser_t = typename precision_chooser<p>::type;
		
}