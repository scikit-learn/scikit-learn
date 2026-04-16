//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_LANCZOS
#define BOOST_MATH_SPECIAL_FUNCTIONS_LANCZOS

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/config.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/mpl/if.hpp>
#include <boost/limits.hpp>
#include <boost/cstdint.hpp>
#include <boost/math/tools/rational.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/mpl/less_equal.hpp>

#include <limits.h>

#if defined(__GNUC__) && defined(BOOST_MATH_USE_FLOAT128)
//
// This is the only way we can avoid
// warning: non-standard suffix on floating constant [-Wpedantic]
// when building with -Wall -pedantic.  Neither __extension__
// nor #pragma dianostic ignored work :(
//
#pragma GCC system_header
#endif

namespace boost{ namespace math{ namespace lanczos{

//
// Individual lanczos approximations start here.
//
// Optimal values for G for each N are taken from
// http://web.mala.bc.ca/pughg/phdThesis/phdThesis.pdf,
// as are the theoretical error bounds.
//
// Constants calculated using the method described by Godfrey
// http://my.fit.edu/~gabdo/gamma.txt and elaborated by Toth at
// http://www.rskey.org/gamma.htm using NTL::RR at 1000 bit precision.
//
// Begin with a small helper to force initialization of constants prior
// to main.  This makes the constant initialization thread safe, even
// when called with a user-defined number type.
//
template <class Lanczos, class T>
struct lanczos_initializer
{
   struct init
   {
      init()
      {
         T t(1);
         Lanczos::lanczos_sum(t);
         Lanczos::lanczos_sum_expG_scaled(t);
         Lanczos::lanczos_sum_near_1(t);
         Lanczos::lanczos_sum_near_2(t);
         Lanczos::g();
      }
      void force_instantiate()const{}
   };
   static const init initializer;
   static void force_instantiate()
   {
      initializer.force_instantiate();
   }
};
template <class Lanczos, class T>
typename lanczos_initializer<Lanczos, T>::init const lanczos_initializer<Lanczos, T>::initializer;
//
// Lanczos Coefficients for N=6 G=5.581
// Max experimental error (with arbitary precision arithmetic) 9.516e-12
// Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006
//
struct lanczos6 : public mpl::int_<35>
{
   //
   // Produces slightly better than float precision when evaluated at
   // double precision:
   //
   template <class T>
   static T lanczos_sum(const T& z)
   {
      lanczos_initializer<lanczos6, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[6] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 8706.349592549009182288174442774377925882)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 8523.650341121874633477483696775067709735)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 3338.029219476423550899999750161289306564)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 653.6424994294008795995653541449610986791)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 63.99951844938187085666201263218840287667)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 2.506628274631006311133031631822390264407))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint16_t) denom[6] = {
         static_cast<boost::uint16_t>(0u),
         static_cast<boost::uint16_t>(24u),
         static_cast<boost::uint16_t>(50u),
         static_cast<boost::uint16_t>(35u),
         static_cast<boost::uint16_t>(10u),
         static_cast<boost::uint16_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      lanczos_initializer<lanczos6, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[6] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 32.81244541029783471623665933780748627823)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 32.12388941444332003446077108933558534361)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 12.58034729455216106950851080138931470954)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 2.463444478353241423633780693218408889251)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 0.2412010548258800231126240760264822486599)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 0.009446967704539249494420221613134244048319))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint16_t) denom[6] = {
         static_cast<boost::uint16_t>(0u),
         static_cast<boost::uint16_t>(24u),
         static_cast<boost::uint16_t>(50u),
         static_cast<boost::uint16_t>(35u),
         static_cast<boost::uint16_t>(10u),
         static_cast<boost::uint16_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      lanczos_initializer<lanczos6, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[5] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 2.044879010930422922760429926121241330235)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, -2.751366405578505366591317846728753993668)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 1.02282965224225004296750609604264824677)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, -0.09786124911582813985028889636665335893627)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 0.0009829742267506615183144364420540766510112)),
      };
      T result = 0;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(k*dz + k*k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      lanczos_initializer<lanczos6, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[5] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 5.748142489536043490764289256167080091892)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, -7.734074268282457156081021756682138251825)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 2.875167944990511006997713242805893543947)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, -0.2750873773533504542306766137703788781776)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 35, 0.002763134585812698552178368447708846850353)),
      };
      T result = 0;
      T z = dz + 2;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(z + k*z + k*k - 1);
      }
      return result;
   }

   static double g(){ return 5.581000000000000405009359383257105946541; }
};

//
// Lanczos Coefficients for N=11 G=10.900511
// Max experimental error (with arbitary precision arithmetic) 2.16676e-19
// Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006
//
struct lanczos11 : public mpl::int_<60>
{
   //
   // Produces slightly better than double precision when evaluated at
   // extended-double precision:
   //
   template <class T>
   static T lanczos_sum(const T& z)
   {
      lanczos_initializer<lanczos11, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[11] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 38474670393.31776828316099004518914832218)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 36857665043.51950660081971227404959150474)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 15889202453.72942008945006665994637853242)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 4059208354.298834770194507810788393801607)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 680547661.1834733286087695557084801366446)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 78239755.00312005289816041245285376206263)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 6246580.776401795264013335510453568106366)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 341986.3488721347032223777872763188768288)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 12287.19451182455120096222044424100527629)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 261.6140441641668190791708576058805625502)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 2.506628274631000502415573855452633787834))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint32_t) denom[11] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(362880u),
         static_cast<boost::uint32_t>(1026576u),
         static_cast<boost::uint32_t>(1172700u),
         static_cast<boost::uint32_t>(723680u),
         static_cast<boost::uint32_t>(269325u),
         static_cast<boost::uint32_t>(63273u),
         static_cast<boost::uint32_t>(9450u),
         static_cast<boost::uint32_t>(870u),
         static_cast<boost::uint32_t>(45u),
         static_cast<boost::uint32_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      lanczos_initializer<lanczos11, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[11] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 709811.662581657956893540610814842699825)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 679979.847415722640161734319823103390728)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 293136.785721159725251629480984140341656)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 74887.5403291467179935942448101441897121)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 12555.29058241386295096255111537516768137)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 1443.42992444170669746078056942194198252)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 115.2419459613734722083208906727972935065)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 6.30923920573262762719523981992008976989)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 0.2266840463022436475495508977579735223818)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 0.004826466289237661857584712046231435101741)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 0.4624429436045378766270459638520555557321e-4))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint32_t) denom[11] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(362880u),
         static_cast<boost::uint32_t>(1026576u),
         static_cast<boost::uint32_t>(1172700u),
         static_cast<boost::uint32_t>(723680u),
         static_cast<boost::uint32_t>(269325u),
         static_cast<boost::uint32_t>(63273u),
         static_cast<boost::uint32_t>(9450u),
         static_cast<boost::uint32_t>(870u),
         static_cast<boost::uint32_t>(45u),
         static_cast<boost::uint32_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      lanczos_initializer<lanczos11, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[10] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 4.005853070677940377969080796551266387954)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -13.17044315127646469834125159673527183164)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 17.19146865350790353683895137079288129318)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -11.36446409067666626185701599196274701126)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 4.024801119349323770107694133829772634737)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -0.7445703262078094128346501724255463005006)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 0.06513861351917497265045550019547857713172)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -0.00217899958561830354633560009312512312758)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 0.17655204574495137651670832229571934738e-4)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -0.1036282091079938047775645941885460820853e-7)),
      };
      T result = 0;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(k*dz + k*k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      lanczos_initializer<lanczos11, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[10] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 19.05889633808148715159575716844556056056)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -62.66183664701721716960978577959655644762)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 81.7929198065004751699057192860287512027)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -54.06941772964234828416072865069196553015)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 19.14904664790693019642068229478769661515)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -3.542488556926667589704590409095331790317)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 0.3099140334815639910894627700232804503017)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -0.01036716187296241640634252431913030440825)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, 0.8399926504443119927673843789048514017761e-4)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 60, -0.493038376656195010308610694048822561263e-7)),
      };
      T result = 0;
      T z = dz + 2;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(z + k*z + k*k - 1);
      }
      return result;
   }

   static double g(){ return 10.90051099999999983936049829935654997826; }
};

//
// Lanczos Coefficients for N=13 G=13.144565
// Max experimental error (with arbitary precision arithmetic) 9.2213e-23
// Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006
//
struct lanczos13 : public mpl::int_<72>
{
   //
   // Produces slightly better than extended-double precision when evaluated at
   // higher precision:
   //
   template <class T>
   static T lanczos_sum(const T& z)
   {
      lanczos_initializer<lanczos13, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[13] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 44012138428004.60895436261759919070125699)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 41590453358593.20051581730723108131357995)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 18013842787117.99677796276038389462742949)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 4728736263475.388896889723995205703970787)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 837910083628.4046470415724300225777912264)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 105583707273.4299344907359855510105321192)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 9701363618.494999493386608345339104922694)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 654914397.5482052641016767125048538245644)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 32238322.94213356530668889463945849409184)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 1128514.219497091438040721811544858643121)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 26665.79378459858944762533958798805525125)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 381.8801248632926870394389468349331394196)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 2.506628274631000502415763426076722427007))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint32_t) denom[13] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(39916800u),
         static_cast<boost::uint32_t>(120543840u),
         static_cast<boost::uint32_t>(150917976u),
         static_cast<boost::uint32_t>(105258076u),
         static_cast<boost::uint32_t>(45995730u),
         static_cast<boost::uint32_t>(13339535u),
         static_cast<boost::uint32_t>(2637558u),
         static_cast<boost::uint32_t>(357423u),
         static_cast<boost::uint32_t>(32670u),
         static_cast<boost::uint32_t>(1925u),
         static_cast<boost::uint32_t>(66u),
         static_cast<boost::uint32_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      lanczos_initializer<lanczos13, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[13] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 86091529.53418537217994842267760536134841)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 81354505.17858011242874285785316135398567)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 35236626.38815461910817650960734605416521)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 9249814.988024471294683815872977672237195)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 1639024.216687146960253839656643518985826)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 206530.8157641225032631778026076868855623)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 18976.70193530288915698282139308582105936)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 1281.068909912559479885759622791374106059)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 63.06093343420234536146194868906771599354)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 2.207470909792527638222674678171050209691)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 0.05216058694613505427476207805814960742102)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 0.0007469903808915448316510079585999893674101)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 0.4903180573459871862552197089738373164184e-5))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint32_t) denom[13] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(39916800u),
         static_cast<boost::uint32_t>(120543840u),
         static_cast<boost::uint32_t>(150917976u),
         static_cast<boost::uint32_t>(105258076u),
         static_cast<boost::uint32_t>(45995730u),
         static_cast<boost::uint32_t>(13339535u),
         static_cast<boost::uint32_t>(2637558u),
         static_cast<boost::uint32_t>(357423u),
         static_cast<boost::uint32_t>(32670u),
         static_cast<boost::uint32_t>(1925u),
         static_cast<boost::uint32_t>(66u),
         static_cast<boost::uint32_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      lanczos_initializer<lanczos13, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[12] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 4.832115561461656947793029596285626840312)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -19.86441536140337740383120735104359034688)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 33.9927422807443239927197864963170585331)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -31.41520692249765980987427413991250886138)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 17.0270866009599345679868972409543597821)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -5.5077216950865501362506920516723682167)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 1.037811741948214855286817963800439373362)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -0.106640468537356182313660880481398642811)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 0.005276450526660653288757565778182586742831)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -0.0001000935625597121545867453746252064770029)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 0.462590910138598083940803704521211569234e-6)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -0.1735307814426389420248044907765671743012e-9)),
      };
      T result = 0;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(k*dz + k*k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      lanczos_initializer<lanczos13, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[12] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 26.96979819614830698367887026728396466395)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -110.8705424709385114023884328797900204863)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 189.7258846119231466417015694690434770085)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -175.3397202971107486383321670769397356553)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 95.03437648691551457087250340903980824948)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -30.7406022781665264273675797983497141978)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 5.792405601630517993355102578874590410552)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -0.5951993240669148697377539518639997795831)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 0.02944979359164017509944724739946255067671)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -0.0005586586555377030921194246330399163602684)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, 0.2581888478270733025288922038673392636029e-5)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 72, -0.9685385411006641478305219367315965391289e-9)),
      };
      T result = 0;
      T z = dz + 2;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(z + k*z + k*k - 1);
      }
      return result;
   }

   static double g(){ return 13.1445650000000000545696821063756942749; }
};

//
// Lanczos Coefficients for N=22 G=22.61891
// Max experimental error (with arbitary precision arithmetic) 2.9524e-38
// Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006
//
struct lanczos22 : public mpl::int_<120>
{
   //
   // Produces slightly better than 128-bit long-double precision when 
   // evaluated at higher precision:
   //
   template <class T>
   static T lanczos_sum(const T& z)
   {
      lanczos_initializer<lanczos22, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[22] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 46198410803245094237463011094.12173081986)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 43735859291852324413622037436.321513777)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 19716607234435171720534556386.97481377748)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 5629401471315018442177955161.245623932129)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 1142024910634417138386281569.245580222392)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 175048529315951173131586747.695329230778)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 21044290245653709191654675.41581372963167)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 2033001410561031998451380.335553678782601)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 160394318862140953773928.8736211601848891)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 10444944438396359705707.48957290388740896)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 565075825801617290121.1466393747967538948)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 25475874292116227538.99448534450411942597)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 957135055846602154.6720835535232270205725)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 29874506304047462.23662392445173880821515)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 769651310384737.2749087590725764959689181)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 16193289100889.15989633624378404096011797)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 273781151680.6807433264462376754578933261)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 3630485900.32917021712188739762161583295)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 36374352.05577334277856865691538582936484)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 258945.7742115532455441786924971194951043)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 1167.501919472435718934219997431551246996)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 2.50662827463100050241576528481104525333))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint64_t) denom[22] = {
         BOOST_MATH_INT_VALUE_SUFFIX(0, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(2432902008176640000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(8752948036761600000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(13803759753640704000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(12870931245150988800, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(8037811822645051776, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(3599979517947607200, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1206647803780373360, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(311333643161390640, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(63030812099294896, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(10142299865511450, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1307535010540395, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(135585182899530, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(11310276995381, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(756111184500, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(40171771630, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1672280820, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(53327946, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1256850, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(20615, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(210, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1, uLL)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      lanczos_initializer<lanczos22, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[22] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 6939996264376682180.277485395074954356211)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 6570067992110214451.87201438870245659384)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 2961859037444440551.986724631496417064121)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 845657339772791245.3541226499766163431651)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 171556737035449095.2475716923888737881837)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 26296059072490867.7822441885603400926007)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 3161305619652108.433798300149816829198706)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 305400596026022.4774396904484542582526472)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 24094681058862.55120507202622377623528108)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 1569055604375.919477574824168939428328839)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 84886558909.02047889339710230696942513159)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 3827024985.166751989686050643579753162298)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 143782298.9273215199098728674282885500522)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 4487794.24541641841336786238909171265944)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 115618.2025760830513505888216285273541959)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 2432.580773108508276957461757328744780439)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 41.12782532742893597168530008461874360191)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.5453771709477689805460179187388702295792)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.005464211062612080347167337964166505282809)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.388992321263586767037090706042788910953e-4)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.1753839324538447655939518484052327068859e-6)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.3765495513732730583386223384116545391759e-9))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint64_t) denom[22] = {
         BOOST_MATH_INT_VALUE_SUFFIX(0, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(2432902008176640000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(8752948036761600000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(13803759753640704000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(12870931245150988800, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(8037811822645051776, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(3599979517947607200, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1206647803780373360, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(311333643161390640, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(63030812099294896, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(10142299865511450, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1307535010540395, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(135585182899530, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(11310276995381, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(756111184500, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(40171771630, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1672280820, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(53327946, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1256850, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(20615, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(210, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1, uLL)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      lanczos_initializer<lanczos22, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[21] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 8.318998691953337183034781139546384476554)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -63.15415991415959158214140353299240638675)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 217.3108224383632868591462242669081540163)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -448.5134281386108366899784093610397354889)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 619.2903759363285456927248474593012711346)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -604.1630177420625418522025080080444177046)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 428.8166750424646119935047118287362193314)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -224.6988753721310913866347429589434550302)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 87.32181627555510833499451817622786940961)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -25.07866854821128965662498003029199058098)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 5.264398125689025351448861011657789005392)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -0.792518936256495243383586076579921559914)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.08317448364744713773350272460937904691566)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -0.005845345166274053157781068150827567998882)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.0002599412126352082483326238522490030412391)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -0.6748102079670763884917431338234783496303e-5)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.908824383434109002762325095643458603605e-7)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -0.5299325929309389890892469299969669579725e-9)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.994306085859549890267983602248532869362e-12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -0.3499893692975262747371544905820891835298e-15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.7260746353663365145454867069182884694961e-20)),
      };
      T result = 0;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(k*dz + k*k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      lanczos_initializer<lanczos22, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[21] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 75.39272007105208086018421070699575462226)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -572.3481967049935412452681346759966390319)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 1969.426202741555335078065370698955484358)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -4064.74968778032030891520063865996757519)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 5612.452614138013929794736248384309574814)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -5475.357667500026172903620177988213902339)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 3886.243614216111328329547926490398103492)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -2036.382026072125407192448069428134470564)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 791.3727954936062108045551843636692287652)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -227.2808432388436552794021219198885223122)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 47.70974355562144229897637024320739257284)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -7.182373807798293545187073539819697141572)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.7537866989631514559601547530490976100468)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -0.05297470142240154822658739758236594717787)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.00235577330936380542539812701472320434133)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -0.6115613067659273118098229498679502138802e-4)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.8236417010170941915758315020695551724181e-6)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -0.4802628430993048190311242611330072198089e-8)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.9011113376981524418952720279739624707342e-11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, -0.3171854152689711198382455703658589996796e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 120, 0.6580207998808093935798753964580596673177e-19)),
      };
      T result = 0;
      T z = dz + 2;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(z + k*z + k*k - 1);
      }
      return result;
   }

   static double g(){ return 22.61890999999999962710717227309942245483; }
};

//
// Lanczos Coefficients for N=6 G=1.428456135094165802001953125
// Max experimental error (with arbitary precision arithmetic) 8.111667e-8
// Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006
//
struct lanczos6m24 : public mpl::int_<24>
{
   //
   // Use for float precision, when evaluated as a float:
   //
   template <class T>
   static T lanczos_sum(const T& z)
   {
      static const T num[6] = {
         static_cast<T>(58.52061591769095910314047740215847630266L),
         static_cast<T>(182.5248962595894264831189414768236280862L),
         static_cast<T>(211.0971093028510041839168287718170827259L),
         static_cast<T>(112.2526547883668146736465390902227161763L),
         static_cast<T>(27.5192015197455403062503721613097825345L),
         static_cast<T>(2.50662858515256974113978724717473206342L)
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint16_t) denom[6] = {
         static_cast<boost::uint16_t>(0u),
         static_cast<boost::uint16_t>(24u),
         static_cast<boost::uint16_t>(50u),
         static_cast<boost::uint16_t>(35u),
         static_cast<boost::uint16_t>(10u),
         static_cast<boost::uint16_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      static const T num[6] = {
         static_cast<T>(14.0261432874996476619570577285003839357L),
         static_cast<T>(43.74732405540314316089531289293124360129L),
         static_cast<T>(50.59547402616588964511581430025589038612L),
         static_cast<T>(26.90456680562548195593733429204228910299L),
         static_cast<T>(6.595765571169314946316366571954421695196L),
         static_cast<T>(0.6007854010515290065101128585795542383721L)
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint16_t) denom[6] = {
         static_cast<boost::uint16_t>(0u),
         static_cast<boost::uint16_t>(24u),
         static_cast<boost::uint16_t>(50u),
         static_cast<boost::uint16_t>(35u),
         static_cast<boost::uint16_t>(10u),
         static_cast<boost::uint16_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      static const T d[5] = {
         static_cast<T>(0.4922488055204602807654354732674868442106L),
         static_cast<T>(0.004954497451132152436631238060933905650346L),
         static_cast<T>(-0.003374784572167105840686977985330859371848L),
         static_cast<T>(0.001924276018962061937026396537786414831385L),
         static_cast<T>(-0.00056533046336427583708166383712907694434L),
      };
      T result = 0;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(k*dz + k*k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      static const T d[5] = {
         static_cast<T>(0.6534966888520080645505805298901130485464L),
         static_cast<T>(0.006577461728560758362509168026049182707101L),
         static_cast<T>(-0.004480276069269967207178373559014835978161L),
         static_cast<T>(0.00255461870648818292376982818026706528842L),
         static_cast<T>(-0.000750517993690428370380996157470900204524L),
      };
      T result = 0;
      T z = dz + 2;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(z + k*z + k*k - 1);
      }
      return result;
   }

   static double g(){ return 1.428456135094165802001953125; }
};

//
// Lanczos Coefficients for N=13 G=6.024680040776729583740234375
// Max experimental error (with arbitary precision arithmetic) 1.196214e-17
// Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006
//
struct lanczos13m53 : public mpl::int_<53>
{
   //
   // Use for double precision, when evaluated as a double:
   //
   template <class T>
   static T lanczos_sum(const T& z)
   {
      static const T num[13] = {
         static_cast<T>(23531376880.41075968857200767445163675473L),
         static_cast<T>(42919803642.64909876895789904700198885093L),
         static_cast<T>(35711959237.35566804944018545154716670596L),
         static_cast<T>(17921034426.03720969991975575445893111267L),
         static_cast<T>(6039542586.35202800506429164430729792107L),
         static_cast<T>(1439720407.311721673663223072794912393972L),
         static_cast<T>(248874557.8620541565114603864132294232163L),
         static_cast<T>(31426415.58540019438061423162831820536287L),
         static_cast<T>(2876370.628935372441225409051620849613599L),
         static_cast<T>(186056.2653952234950402949897160456992822L),
         static_cast<T>(8071.672002365816210638002902272250613822L),
         static_cast<T>(210.8242777515793458725097339207133627117L),
         static_cast<T>(2.506628274631000270164908177133837338626L)
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint32_t) denom[13] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(39916800u),
         static_cast<boost::uint32_t>(120543840u),
         static_cast<boost::uint32_t>(150917976u),
         static_cast<boost::uint32_t>(105258076u),
         static_cast<boost::uint32_t>(45995730u),
         static_cast<boost::uint32_t>(13339535u),
         static_cast<boost::uint32_t>(2637558u),
         static_cast<boost::uint32_t>(357423u),
         static_cast<boost::uint32_t>(32670u),
         static_cast<boost::uint32_t>(1925u),
         static_cast<boost::uint32_t>(66u),
         static_cast<boost::uint32_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      static const T num[13] = {
         static_cast<T>(56906521.91347156388090791033559122686859L),
         static_cast<T>(103794043.1163445451906271053616070238554L),
         static_cast<T>(86363131.28813859145546927288977868422342L),
         static_cast<T>(43338889.32467613834773723740590533316085L),
         static_cast<T>(14605578.08768506808414169982791359218571L),
         static_cast<T>(3481712.15498064590882071018964774556468L),
         static_cast<T>(601859.6171681098786670226533699352302507L),
         static_cast<T>(75999.29304014542649875303443598909137092L),
         static_cast<T>(6955.999602515376140356310115515198987526L),
         static_cast<T>(449.9445569063168119446858607650988409623L),
         static_cast<T>(19.51992788247617482847860966235652136208L),
         static_cast<T>(0.5098416655656676188125178644804694509993L),
         static_cast<T>(0.006061842346248906525783753964555936883222L)
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint32_t) denom[13] = {
         static_cast<boost::uint32_t>(0u),
         static_cast<boost::uint32_t>(39916800u),
         static_cast<boost::uint32_t>(120543840u),
         static_cast<boost::uint32_t>(150917976u),
         static_cast<boost::uint32_t>(105258076u),
         static_cast<boost::uint32_t>(45995730u),
         static_cast<boost::uint32_t>(13339535u),
         static_cast<boost::uint32_t>(2637558u),
         static_cast<boost::uint32_t>(357423u),
         static_cast<boost::uint32_t>(32670u),
         static_cast<boost::uint32_t>(1925u),
         static_cast<boost::uint32_t>(66u),
         static_cast<boost::uint32_t>(1u)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      static const T d[12] = {
         static_cast<T>(2.208709979316623790862569924861841433016L),
         static_cast<T>(-3.327150580651624233553677113928873034916L),
         static_cast<T>(1.483082862367253753040442933770164111678L),
         static_cast<T>(-0.1993758927614728757314233026257810172008L),
         static_cast<T>(0.004785200610085071473880915854204301886437L),
         static_cast<T>(-0.1515973019871092388943437623825208095123e-5L),
         static_cast<T>(-0.2752907702903126466004207345038327818713e-7L),
         static_cast<T>(0.3075580174791348492737947340039992829546e-7L),
         static_cast<T>(-0.1933117898880828348692541394841204288047e-7L),
         static_cast<T>(0.8690926181038057039526127422002498960172e-8L),
         static_cast<T>(-0.2499505151487868335680273909354071938387e-8L),
         static_cast<T>(0.3394643171893132535170101292240837927725e-9L),
      };
      T result = 0;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(k*dz + k*k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      static const T d[12] = {
         static_cast<T>(6.565936202082889535528455955485877361223L),
         static_cast<T>(-9.8907772644920670589288081640128194231L),
         static_cast<T>(4.408830289125943377923077727900630927902L),
         static_cast<T>(-0.5926941084905061794445733628891024027949L),
         static_cast<T>(0.01422519127192419234315002746252160965831L),
         static_cast<T>(-0.4506604409707170077136555010018549819192e-5L),
         static_cast<T>(-0.8183698410724358930823737982119474130069e-7L),
         static_cast<T>(0.9142922068165324132060550591210267992072e-7L),
         static_cast<T>(-0.5746670642147041587497159649318454348117e-7L),
         static_cast<T>(0.2583592566524439230844378948704262291927e-7L),
         static_cast<T>(-0.7430396708998719707642735577238449585822e-8L),
         static_cast<T>(0.1009141566987569892221439918230042368112e-8L),
      };
      T result = 0;
      T z = dz + 2;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(z + k*z + k*k - 1);
      }
      return result;
   }

   static double g(){ return 6.024680040776729583740234375; }
};

//
// Lanczos Coefficients for N=17 G=12.2252227365970611572265625
// Max experimental error (with arbitary precision arithmetic) 2.7699e-26
// Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006
//
struct lanczos17m64 : public mpl::int_<64>
{
   //
   // Use for extended-double precision, when evaluated as an extended-double:
   //
   template <class T>
   static T lanczos_sum(const T& z)
   {
      lanczos_initializer<lanczos17m64, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[17] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 553681095419291969.2230556393350368550504)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 731918863887667017.2511276782146694632234)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 453393234285807339.4627124634539085143364)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 174701893724452790.3546219631779712198035)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 46866125995234723.82897281620357050883077)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 9281280675933215.169109622777099699054272)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 1403600894156674.551057997617468721789536)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 165345984157572.7305349809894046783973837)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 15333629842677.31531822808737907246817024)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 1123152927963.956626161137169462874517318)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 64763127437.92329018717775593533620578237)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 2908830362.657527782848828237106640944457)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 99764700.56999856729959383751710026787811)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 2525791.604886139959837791244686290089331)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 44516.94034970167828580039370201346554872)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 488.0063567520005730476791712814838113252)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 2.50662827463100050241576877135758834683))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint64_t) denom[17] = {
         BOOST_MATH_INT_VALUE_SUFFIX(0, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1307674368000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(4339163001600, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(6165817614720, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(5056995703824, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(2706813345600, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1009672107080, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(272803210680, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(54631129553, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(8207628000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(928095740, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(78558480, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(4899622, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(218400, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(6580, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(120, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1, uLL)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      lanczos_initializer<lanczos17m64, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[17] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 2715894658327.717377557655133124376674911)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 3590179526097.912105038525528721129550434)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 2223966599737.814969312127353235818710172)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 856940834518.9562481809925866825485883417)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 229885871668.749072933597446453399395469)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 45526171687.54610815813502794395753410032)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 6884887713.165178784550917647709216424823)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 811048596.1407531864760282453852372777439)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 75213915.96540822314499613623119501704812)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 5509245.417224265151697527957954952830126)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 317673.5368435419126714931842182369574221)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 14268.27989845035520147014373320337523596)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 489.3618720403263670213909083601787814792)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 12.38941330038454449295883217865458609584)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.2183627389504614963941574507281683147897)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.002393749522058449186690627996063983095463)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.1229541408909435212800785616808830746135e-4))
      };
      static const BOOST_MATH_INT_TABLE_TYPE(T, boost::uint64_t) denom[17] = {
         BOOST_MATH_INT_VALUE_SUFFIX(0, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1307674368000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(4339163001600, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(6165817614720, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(5056995703824, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(2706813345600, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1009672107080, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(272803210680, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(54631129553, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(8207628000, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(928095740, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(78558480, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(4899622, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(218400, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(6580, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(120, uLL),
         BOOST_MATH_INT_VALUE_SUFFIX(1, uLL)
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      lanczos_initializer<lanczos17m64, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[16] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 4.493645054286536365763334986866616581265)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -16.95716370392468543800733966378143997694)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 26.19196892983737527836811770970479846644)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -21.3659076437988814488356323758179283908)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 9.913992596774556590710751047594507535764)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -2.62888300018780199210536267080940382158)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.3807056693542503606384861890663080735588)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.02714647489697685807340312061034730486958)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.0007815484715461206757220527133967191796747)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.6108630817371501052576880554048972272435e-5)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.5037380238864836824167713635482801545086e-8)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.1483232144262638814568926925964858237006e-13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.1346609158752142460943888149156716841693e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.660492688923978805315914918995410340796e-15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.1472114697343266749193617793755763792681e-15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.1410901942033374651613542904678399264447e-16)),
      };
      T result = 0;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(k*dz + k*k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      lanczos_initializer<lanczos17m64, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[16] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 23.56409085052261327114594781581930373708)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -88.92116338946308797946237246006238652361)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 137.3472822086847596961177383569603988797)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -112.0400438263562152489272966461114852861)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 51.98768915202973863076166956576777843805)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -13.78552090862799358221343319574970124948)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 1.996371068830872830250406773917646121742)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.1423525874909934506274738563671862576161)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.004098338646046865122459664947239111298524)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.3203286637326511000882086573060433529094e-4)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.2641536751640138646146395939004587594407e-7)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.7777876663062235617693516558976641009819e-13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.7061443477097101636871806229515157914789e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.3463537849537988455590834887691613484813e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, 0.7719578215795234036320348283011129450595e-15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 64, -0.7398586479708476329563577384044188912075e-16)),
      };
      T result = 0;
      T z = dz + 2;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(z + k*z + k*k - 1);
      }
      return result;
   }

   static double g(){ return 12.2252227365970611572265625; }
};

//
// Lanczos Coefficients for N=24 G=20.3209821879863739013671875
// Max experimental error (with arbitary precision arithmetic) 1.0541e-38
// Generated with compiler: Microsoft Visual C++ version 8.0 on Win32 at Mar 23 2006
//
struct lanczos24m113 : public mpl::int_<113>
{
   //
   // Use for long-double precision, when evaluated as an long-double:
   //
   template <class T>
   static T lanczos_sum(const T& z)
   {
      lanczos_initializer<lanczos24m113, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[24] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 2029889364934367661624137213253.22102954656825019111612712252027267955023987678816620961507)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 2338599599286656537526273232565.2727349714338768161421882478417543004440597874814359063158)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1288527989493833400335117708406.3953711906175960449186720680201425446299360322830739180195)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 451779745834728745064649902914.550539158066332484594436145043388809847364393288132164411521)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 113141284461097964029239556815.291212318665536114012605167994061291631013303788706545334708)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 21533689802794625866812941616.7509064680880468667055339259146063256555368135236149614592432)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 3235510315314840089932120340.71494940111731241353655381919722177496659303550321056514776757)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 393537392344185475704891959.081297108513472083749083165179784098220158201055270548272414314)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 39418265082950435024868801.5005452240816902251477336582325944930252142622315101857742955673)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 3290158764187118871697791.05850632319194734270969161036889516414516566453884272345518372696)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 230677110449632078321772.618245845856640677845629174549731890660612368500786684333975350954)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 13652233645509183190158.5916189185218250859402806777406323001463296297553612462737044693697)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 683661466754325350495.216655026531202476397782296585200982429378069417193575896602446904762)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 28967871782219334117.0122379171041074970463982134039409352925258212207710168851968215545064)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1036104088560167006.2022834098572346459442601718514554488352117620272232373622553429728555)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 31128490785613152.8380102669349814751268126141105475287632676569913936040772990253369753962)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 779327504127342.536207878988196814811198475410572992436243686674896894543126229424358472541)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 16067543181294.643350688789124777020407337133926174150582333950666044399234540521336771876)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 268161795520.300916569439413185778557212729611517883948634711190170998896514639936969855484)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 3533216359.10528191668842486732408440112703691790824611391987708562111396961696753452085068)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 35378979.5479656110614685178752543826919239614088343789329169535932709470588426584501652577)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 253034.881362204346444503097491737872930637147096453940375713745904094735506180552724766444)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1151.61895453463992438325318456328526085882924197763140514450975619271382783957699017875304)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 2.50662827463100050241576528481104515966515623051532908941425544355490413900497467936202516))
      };
      static const T denom[24] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.112400072777760768e22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.414847677933545472e22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 6756146673770930688000.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 6548684852703068697600.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 4280722865357147142912.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 2021687376910682741568.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 720308216440924653696.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 199321978221066137360.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 43714229649594412832.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 7707401101297361068.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1103230881185949736.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 129006659818331295.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 12363045847086207.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 971250460939913.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 62382416421941.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 3256091103430.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 136717357942.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 4546047198.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 116896626.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 2240315.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 30107.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 253.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1.0))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      lanczos_initializer<lanczos24m113, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T num[24] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 3035162425359883494754.02878223286972654682199012688209026810841953293372712802258398358538)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 3496756894406430103600.16057175075063458536101374170860226963245118484234495645518505519827)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1926652656689320888654.01954015145958293168365236755537645929361841917596501251362171653478)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 675517066488272766316.083023742440619929434602223726894748181327187670231286180156444871912)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 169172853104918752780.086262749564831660238912144573032141700464995906149421555926000038492)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 32197935167225605785.6444116302160245528783954573163541751756353183343357329404208062043808)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 4837849542714083249.37587447454818124327561966323276633775195138872820542242539845253171632)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 588431038090493242.308438203986649553459461798968819276505178004064031201740043314534404158)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 58939585141634058.6206417889192563007809470547755357240808035714047014324843817783741669733)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 4919561837722192.82991866530802080996138070630296720420704876654726991998309206256077395868)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 344916580244240.407442753122831512004021081677987651622305356145640394384006997569631719101)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 20413302960687.8250598845969238472629322716685686993835561234733641729957841485003560103066)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1022234822943.78400752460970689311934727763870970686747383486600540378889311406851534545789)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 43313787191.9821354846952908076307094286897439975815501673706144217246093900159173598852503)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1549219505.59667418528481770869280437577581951167003505825834192510436144666564648361001914)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 46544421.1998761919380541579358096705925369145324466147390364674998568485110045455014967149)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1165278.06807504975090675074910052763026564833951579556132777702952882101173607903881127542)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 24024.759267256769471083727721827405338569868270177779485912486668586611981795179894572115)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 400.965008113421955824358063769761286758463521789765880962939528760888853281920872064838918)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 5.28299015654478269617039029170846385138134929147421558771949982217659507918482272439717603)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.0528999024412510102409256676599360516359062802002483877724963720047531347449011629466149805)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.000378346710654740685454266569593414561162134092347356968516522170279688139165340746957511115)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.172194142179211139195966608011235161516824700287310869949928393345257114743230967204370963e-5)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.374799931707148855771381263542708435935402853962736029347951399323367765509988401336565436e-8))
      };
      static const T denom[24] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.112400072777760768e22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.414847677933545472e22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 6756146673770930688000.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 6548684852703068697600.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 4280722865357147142912.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 2021687376910682741568.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 720308216440924653696.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 199321978221066137360.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 43714229649594412832.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 7707401101297361068.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1103230881185949736.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 129006659818331295.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 12363045847086207.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 971250460939913.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 62382416421941.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 3256091103430.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 136717357942.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 4546047198.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 116896626.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 2240315.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 30107.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 253.0)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1.0))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      lanczos_initializer<lanczos24m113, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[23] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 7.4734083002469026177867421609938203388868806387315406134072298925733950040583068760685908)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -50.4225805042247530267317342133388132970816607563062253708655085754357843064134941138154171)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 152.288200621747008570784082624444625293884063492396162110698238568311211546361189979357019)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -271.894959539150384169327513139846971255640842175739337449692360299099322742181325023644769)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 319.240102980202312307047586791116902719088581839891008532114107693294261542869734803906793)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -259.493144143048088289689500935518073716201741349569864988870534417890269467336454358361499)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 149.747518319689708813209645403067832020714660918583227716408482877303972685262557460145835)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -61.9261301009341333289187201425188698128684426428003249782448828881580630606817104372760037)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 18.3077524177286961563937379403377462608113523887554047531153187277072451294845795496072365)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -3.82011322251948043097070160584761236869363471824695092089556195047949392738162970152230254)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.549382685505691522516705902336780999493262538301283190963770663549981309645795228539620711)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.0524814679715180697633723771076668718265358076235229045603747927518423453658004287459638024)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.00315392664003333528534120626687784812050217700942910879712808180705014754163256855643360698)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.000110098373127648510519799564665442121339511198561008748083409549601095293123407080388658329)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.19809382866681658224945717689377373458866950897791116315219376038432014207446832310901893e-5)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.152278977408600291408265615203504153130482270424202400677280558181047344681214058227949755e-7)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.364344768076106268872239259083188037615571711218395765792787047015406264051536972018235217e-10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.148897510480440424971521542520683536298361220674662555578951242811522959610991621951203526e-13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.261199241161582662426512749820666625442516059622425213340053324061794752786482115387573582e-18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.780072664167099103420998436901014795601783313858454665485256897090476089641613851903791529e-24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.303465867587106629530056603454807425512962762653755513440561256044986695349304176849392735e-24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.615420597971283870342083342286977366161772327800327789325710571275345878439656918541092056e-25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.499641233843540749369110053005439398774706583601830828776209650445427083113181961630763702e-26)),
      };
      T result = 0;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(k*dz + k*k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      lanczos_initializer<lanczos24m113, T>::force_instantiate(); // Ensure our constants get initialized before main()
      static const T d[23] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 61.4165001061101455341808888883960361969557848005400286332291451422461117307237198559485365)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -414.372973678657049667308134761613915623353625332248315105320470271523320700386200587519147)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1251.50505818554680171298972755376376836161706773644771875668053742215217922228357204561873)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -2234.43389421602399514176336175766511311493214354568097811220122848998413358085613880612158)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 2623.51647746991904821899989145639147785427273427135380151752779100215839537090464785708684)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -2132.51572435428751962745870184529534443305617818870214348386131243463614597272260797772423)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 1230.62572059218405766499842067263311220019173335523810725664442147670956427061920234820189)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -508.90919151163744999377586956023909888833335885805154492270846381061182696305011395981929)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 150.453184562246579758706538566480316921938628645961177699894388251635886834047343195475395)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -31.3937061525822497422230490071156186113405446381476081565548185848237169870395131828731397)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 4.51482916590287954234936829724231512565732528859217337795452389161322923867318809206313688)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.431292919341108177524462194102701868233551186625103849565527515201492276412231365776131952)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.0259189820815586225636729971503340447445001375909094681698918294680345547092233915092128323)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.000904788882557558697594884691337532557729219389814315972435534723829065673966567231504429712)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.162793589759218213439218473348810982422449144393340433592232065020562974405674317564164312e-4)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.125142926178202562426432039899709511761368233479483128438847484617555752948755923647214487e-6)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.299418680048132583204152682950097239197934281178261879500770485862852229898797687301941982e-9)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.122364035267809278675627784883078206654408225276233049012165202996967011873995261617995421e-12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.21465364366598631597052073538883430194257709353929022544344097235100199405814005393447785e-17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.641064035802907518396608051803921688237330857546406669209280666066685733941549058513986818e-23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.249388374622173329690271566855185869111237201309011956145463506483151054813346819490278951e-23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, -0.505752900177513489906064295001851463338022055787536494321532352380960774349054239257683149e-24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 113, 0.410605371184590959139968810080063542546949719163227555918846829816144878123034347778284006e-25)),
      };
      T result = 0;
      T z = dz + 2;
      for(unsigned k = 1; k <= sizeof(d)/sizeof(d[0]); ++k)
      {
         result += (-d[k-1]*dz)/(z + k*z + k*k - 1);
      }
      return result;
   }

   static double g(){ return 20.3209821879863739013671875; }
};

//
// Lanczos Coefficients for N=32 G=2.1471552819013595581054687500000000000000000e+01
// Max experimental error (with 40 digit precision arithmetic) 4.3871855787077623312177313715826599434111453e-39
// Generated with compiler: Microsoft Visual C++ version 14.2 on Win32 at Oct 14 2019
// Type precision was 134 bits or 43 max_digits10
//
struct lanczos32MP : public mpl::int_<134>
{
   template <class T>
   static T lanczos_sum(const T& z)
   {
      static const T num[32] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.5570588792269726580426965328821299006241260e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.1774633762992816326945208100096848962840570e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4729248542870148408945035000721777267825674e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 6.4187209355291098195390045608656819929546438e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.0248306413899929430368029054268103901281039e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 4.9257678361460200143521812997589499752976231e+40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 9.6124719365239930504906027891035869583753288e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.5454747713027643957043517109001371929077089e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.0864696676761302003121425399683010735367223e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.3986590841584046388095075457448938071498037e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.3730272176307717941787236450736722205208718e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.0364763054036420066290166982737191399230330e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.5251577656813796428428395564637209124879141e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.0012588845037002389721073745460292782912025e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 5.7802279869271138404361989987596881149322644e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.9402720678039082125732480764868685957396877e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.3191712125257631604053220568076491847333218e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 5.2201606453728489193256459039591784247070127e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.8201383870559727838136908091985005180163020e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 5.5806105013691586921928194185346247254008596e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4998261429327606587262540663611193050275961e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.5175518216952441498193267690845732518664615e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.1558306683405617426954019692046993902758157e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.2526529449739780296036175515849964203851625e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.8671408513624251160128760408315369701624440e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.3367482193223363697833634493247948824157616e+14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.4092037586229976036163287167214501623150351e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.9923609136198492516415716991452463530279321e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.2704580842864713292739934791084236962038214e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 5.8637257952414489194330242306329389356743799e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.7432900891411285513294395078026439443163899e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.5066282746310005024157652848110452530069853e+00))
      };
      static const T denom[32] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 0.0000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.6525285981219105863630848000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.0596817613895338599493271552000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.9028937852409282099982165606400000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.0707922020245946836608666419200000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.5477949752547197371117812531200000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 8.5189988850542311250318425141248000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.6093078815883681280561453887897600000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.2136536667474513652307465210240000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.3114629767614997850763390570240000000000000e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.4541614716906607001396551577600000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4019376240868075016911422397440000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.2245742324696206305840307600000000000000000e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.0006513636556697864066736800000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.4602661104938986779113940000000000000000000e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.4256361393293766065270064000000000000000000e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.9197210605623737977801375000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.1458832493345014521397750000000000000000000e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.3605580871196332287117500000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.4359416261117272348550000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.4960054586805754087500000000000000000000000e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4090257524223082475000000000000000000000000e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 4.8433867667953267500000000000000000000000000e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4097793282984515000000000000000000000000000e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.4409270792812500000000000000000000000000000e+14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 6.9491892473250000000000000000000000000000000e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.1400943144500000000000000000000000000000000e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4803212690000000000000000000000000000000000e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4631225000000000000000000000000000000000000e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.0338500000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 4.6500000000000000000000000000000000000000000e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.0000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      static const T num[32] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.3676354398462675075959984137356365580662550e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.0303243219777652371982632367487187515581930e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 6.9695330738315109132895792451395673621363046e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.0371873841126513851737915789350097633065208e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 9.5810211111592336765499877083050083223251925e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.3307571834445765987280306477823996158240302e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 4.5483950445870222746510194590139013118104375e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.3128221728465167159919234074383914011057504e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 9.8726824481851158990505892727571391497920715e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.1349889148269349577271512140093856701260605e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.1228605200219708164439121714303738845598879e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 9.6361256470584159095341545227609851733775758e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.2166868932853575755146911619588166896783547e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 4.7377209303689212354445666628772172811909754e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.7350675774071094836814884499810988826239303e+22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.3912674066825924538916027148038869739001663e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 6.2420070983154004111449116807024523016395065e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.4700569185690640515452949188152438288912426e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 8.6124656329989777592191745322776627796111867e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.6406132905056437989704993228158312134178284e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.0968236279962512789710882617045613206867906e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.6644225731453470304909038926424696324288813e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.3859703275819088757835059385783332958674837e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 5.9272583422150418041465884118526583443207265e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 8.8348702102467259218213323557081872404861389e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.1056941535328594779262418139203368431235554e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.1399784061251411774918201133708617903068477e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 9.4273820161750809875685215118516360445304050e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 6.0115080627363704378364976159119115750197259e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.7745767713042344487490571988087795152219347e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 8.2488376091888795863984895909744485000401866e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.1860773896913111295543391748504025973928588e-09))
      };
      static const T denom[32] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 0.0000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.6525285981219105863630848000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.0596817613895338599493271552000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.9028937852409282099982165606400000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.0707922020245946836608666419200000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.5477949752547197371117812531200000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 8.5189988850542311250318425141248000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.6093078815883681280561453887897600000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.2136536667474513652307465210240000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.3114629767614997850763390570240000000000000e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.4541614716906607001396551577600000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4019376240868075016911422397440000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.2245742324696206305840307600000000000000000e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.0006513636556697864066736800000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.4602661104938986779113940000000000000000000e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.4256361393293766065270064000000000000000000e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.9197210605623737977801375000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.1458832493345014521397750000000000000000000e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.3605580871196332287117500000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.4359416261117272348550000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.4960054586805754087500000000000000000000000e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4090257524223082475000000000000000000000000e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 4.8433867667953267500000000000000000000000000e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4097793282984515000000000000000000000000000e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.4409270792812500000000000000000000000000000e+14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 6.9491892473250000000000000000000000000000000e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.1400943144500000000000000000000000000000000e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4803212690000000000000000000000000000000000e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4631225000000000000000000000000000000000000e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.0338500000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 4.6500000000000000000000000000000000000000000e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.0000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      static const T d[31] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 7.8968008940705433227909677660515993038374841e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -5.6618604605116020762015754929520236835404324e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.8292545018854689128055824823953561237213555e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -3.5208766506452851610872857905478685024504482e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 4.4977346993761774005199362646328436255151093e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -4.0215319141694562206816207245107490261482447e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.5868485920551737806269221936902170855620396e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -1.2119342552650938943400026113911660156093453e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 4.1418202444782049081885117917390551996878611e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -1.0248307330153472190143848702750491768975848e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.8061053255438937388879028538451417962537674e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -2.2080708594510418170142862951931773430228843e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.8014009027542132302057253782788244298905648e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -9.2759167304556987616725030218164051987480377e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.7821648880614102293660548885638379600357589e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -4.3159048716588863974468231104830412871370083e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.8802059531699177780993451552740605512217594e-09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -6.0864938064916182822437411537845427325953853e-12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.3077465050230400693183749778781074723249743e-15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -4.4535910089559609770590606433653482492203612e-20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.8463914918024931651359730176750702440539677e-26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -2.6476799036488876225768641543721151210680141e-26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.9859528302828214478942827717619270410965131e-26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -2.4164337794361162374639737594212065825014909e-26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.4264219556008516021575829429243024924344353e-26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -6.3654171220500680155015450276084381767112451e-27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.1771337070979387873887155215691051891368210e-27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -5.6130349210382862255475637378026742528052158e-28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.0357425223049510242618626197816013695138258e-28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -1.2237222565117372586662393874248709023486029e-29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 6.9653890134419774710417169916706146555543512e-31))
      };
      T result = 0;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (k * dz + k * k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      static const T d[31] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 6.8235732906799716623552488416933899152528589e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -4.8923760814221749411260475309979161560716608e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.5806466857096723284428938485933849993255433e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -3.0423661676922237999255010522403805001149858e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.8864627302775745934819737805015890258069096e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -3.4749790611725489418874785394046879079002711e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.2352787155916216728704971718263237763940956e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -1.0472243539148225143721882968062403836164182e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 3.5789194097878442125894349247978066638013227e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -8.8554944097961585263202464603716331582929738e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.5606436359298279012809323203777498460831161e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -1.9079793884371809023191793438331090550554931e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.5565785753912199592703963702155052413697331e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -8.0152581402976896653242934862478024763161144e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.4040502318727807490297102100258877023641629e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -3.7293447818191967510645955265996052825799609e-06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.4887668661450911754046983982188022868005334e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -5.2592989400365993333512441603746877766052195e-11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.9941084528494632503325323388230135305735007e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -3.8483184600922127060811720536956793306527275e-19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.5954546450654998035263505446322855520450213e-25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -2.2878426485811911877035106634321308368820032e-25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 2.5801420414749176901982573044724995471338101e-25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -2.0880244059892646868604082519013410844174147e-25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.2325617535559095736858057307295434296584568e-25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -5.5003147275337403154295953682157801293504851e-26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 1.8812468002260686540475204400230633058447779e-26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -4.8501862565143035234095815929073455137122638e-27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 8.9497824575117543406738759847704395400692055e-28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, -1.0574102876284982539299549303708483236373651e-28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 134, 6.0187464606086614102031360043044718659169199e-30)),
      };
      T result = 0;
      T z = dz + 2;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (z + k * z + k * k - 1);
      }
      return result;
   }

   static double g() { return 2.1471552819013595581054687500000000000000000e+01; }
};

//
// Lanczos Coefficients for N=35 G=2.96640371531248092651367187500000000000000000000000000e+01
// Max experimental error (with 50 digit precision arithmetic) 67eps
// Generated with compiler: Microsoft Visual C++ version 14.2 on Win32 at Oct 14 2019
// Type precision was 168 bits or 53 max_digits10
//
struct lanczos35MP : public mpl::int_<168>
{
   template <class T>
   static T lanczos_sum(const T& z)
   {
      static const T num[35] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.17215050716253100021302249837728942659410271586236104e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.51055117651708470336913962553466820524801246971658127e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.40813458996718289733677017073036013655624930344397267e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.10569518324826607478187974291222641098997506635019681e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.34502197565331471178368569687788687058240547971732391e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.74311603169690571192608960963509140372217014888512918e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 4.50656021978234091874071935392175934984492682009447097e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 6.12703102551730381018400796362603958419580969330315139e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.02844698442195350077632196816248435420923619452768200e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 6.90106767379334717236568166816961185224083190775430842e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.86371531667026447746284883480888667804130713757839681e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 4.34808948517797782155274346690360992144536507118093783e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.83232124439938458545786668616393415008373341980153072e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.62895707563068512468013948922815298700909218398406635e+40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.30384063116420066671650072267242339695473078925159324e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.76258309689585811716178198120267186946262194080905971e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.51837231299916455171135124843484994848995300472356341e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.46324357690180919340289798257560253430931750807924001e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.75333853376321853646128997503611223620394342435525484e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.01719517877315910652307531002686423847077617217874485e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.27861878894319497853745513558138184450369083409359360e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.89640024726662067702004632718605032785787967237099607e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.81537701811791870172286588846619085013138846074815251e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.03090758312551459302562064161308518889144037164899961e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.60538569869661647274451913615710409703905629234367906e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.18176163448730621246454091850022844174919234685832508e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.56586635256765282348264053213197702964352373258511008e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.58289895656990946427745668670352144404744258615044371e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.19373478903102411154024309088124853938046967389531861e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.54192605870424877025476980158698548681325282029269310e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.73027427579217615249706012469272147499107562412573337e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.82675918536460865549992482360500962016208597062710654e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.21869956201943834772161655315196962519434419814106818e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.50897418653428667959996348205296461689142907811767371e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.50662827463100050241576528481104525300698674060984055e+00))
      };
      static const T denom[35] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 0.00000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.68331761881188649551819440128000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.55043336733310191803732770947072000000000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 6.55728779174162547080350866368102400000000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.37352350419052295388404251629977600000000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.72117566475005542296335706764492800000000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.28417720643003773414159612967554252800000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.45822739485943139719482682477713244160000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.16476527817201997988283152951021977600000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.49225481668254064104679479029764121600000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.57726463942545496998486904826347776000000000000000000e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.20859297660335343156864734965859840000000000000000000e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.23364307820330543590375511999050240000000000000000000e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.80750015058176473779293385245398400000000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.28183125026789051815954180232544000000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.49437224233918151570015089338400000000000000000000000e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.37000480501772121324931003824000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.96258640868140652967646352465000000000000000000000000e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.41894262447739018035536664650000000000000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.96452376168568744680811480000000000000000000000000000e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 4.94875410890088264440962800000000000000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.38478815149246067334598000000000000000000000000000000e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.00124085806115519088380000000000000000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.65117470518809938644000000000000000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.15145312544238764840000000000000000000000000000000000e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.12192419709374919000000000000000000000000000000000000e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.22038661704031100000000000000000000000000000000000000e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.40979763670090400000000000000000000000000000000000000e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.29191290647440000000000000000000000000000000000000000e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.04437176604000000000000000000000000000000000000000000e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.21763644400000000000000000000000000000000000000000000e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.60169360000000000000000000000000000000000000000000000e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.51096000000000000000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.61000000000000000000000000000000000000000000000000000e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.00000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      static const T num[35] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.84421398435712762388902267099927585742388886580864424e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.28731583799033736725852757551292030085556435695468295e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.84381150359300352571680869181416248982215282642834936e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 6.68539753215772969226355064737523321566208288321687448e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.76117184320624276162478300964159399462275652881271996e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.59183627116994441494601110756468114877940946273012852e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.90089018057779871758440184258134151304912092733579104e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.02273473587728940068021671629793244969348874651645551e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 9.20304883823127369598764418881022021049206245678741573e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 9.03625836242722113759123056762610636251641913153595812e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.67794913334462808923359541498599600753842936204419932e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.69338859264140114791649895977363900871692586779302150e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.70864158121145435408364940074910197916145829346031858e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.13295647753179115743895667847873122731507276407230715e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.08730493440263847356723847541024859440843056640671533e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 4.92672649809905793239714364398097142490510744815940192e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.98815678372776973689475889094271298156568135487559824e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.15357141696015228406471054927723105303656292491717836e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.29582156512528703674984172534752222415664014582498353e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 6.56951562180494343732211791410530161839249714612303326e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.67422350715677024140556410421772283993277946880053914e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.79254663081905790190270601146772274854974105071798035e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.61465496276608608941993297108655885737613121720232292e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.34987044168298086318822469739196823360923972361455073e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.10209211537761991333937729340544738747931371426736883e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.85679879496413826670691454915567101976631415248412906e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.35974553231926272707704478737590721340254406209650188e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.38204802486455055334129565820015244464343854444712513e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.87247644413155087645140975008088533286977710080244249e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.01899805954981363917258740277358024893572331522514601e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.14314215799519834172753514406176454576793263619287700e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.01075867159821346256470334018168931185179114379271616e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.59576526838074751422330690168945437827562833198707558e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.28525092722679899458094768960179796663588010298597603e-10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.28217919006153582429216342066702743329957749672852350e-13))
      };
      static const T denom[35] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 0.00000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.68331761881188649551819440128000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.55043336733310191803732770947072000000000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 6.55728779174162547080350866368102400000000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.37352350419052295388404251629977600000000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.72117566475005542296335706764492800000000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.28417720643003773414159612967554252800000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.45822739485943139719482682477713244160000000000000000e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.16476527817201997988283152951021977600000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.49225481668254064104679479029764121600000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.57726463942545496998486904826347776000000000000000000e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.20859297660335343156864734965859840000000000000000000e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.23364307820330543590375511999050240000000000000000000e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.80750015058176473779293385245398400000000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.28183125026789051815954180232544000000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.49437224233918151570015089338400000000000000000000000e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.37000480501772121324931003824000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.96258640868140652967646352465000000000000000000000000e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.41894262447739018035536664650000000000000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.96452376168568744680811480000000000000000000000000000e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 4.94875410890088264440962800000000000000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.38478815149246067334598000000000000000000000000000000e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.00124085806115519088380000000000000000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.65117470518809938644000000000000000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.15145312544238764840000000000000000000000000000000000e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.12192419709374919000000000000000000000000000000000000e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.22038661704031100000000000000000000000000000000000000e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.40979763670090400000000000000000000000000000000000000e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.29191290647440000000000000000000000000000000000000000e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.04437176604000000000000000000000000000000000000000000e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.21763644400000000000000000000000000000000000000000000e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.60169360000000000000000000000000000000000000000000000e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.51096000000000000000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.61000000000000000000000000000000000000000000000000000e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.00000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      static const T d[34] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.09112391094335813989230740596164619994797033481760301e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.11095925828443504625574261745581427703630213518975734e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.25795762399049096970854782036824945808150003653799701e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.53596813364794430843820749839089482856771473577902136e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.10150518856225197104044336968686219568570379161186763e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -4.59440015012604275008985164760189204399872934212169290e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.17182126100715726914905437099855137991832224176746144e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -4.52194070233963020921586697582404279809725757428225907e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.11287192546142126817144783716567250328231623970392395e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.70029010284498269942986238607675315541449827223720491e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.39359986548869257712075354365607803246217579350731025e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -2.55857186796151318525331772210231419636656738702520247e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.01976508569781091139272219093293962106866776273944401e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.51622799316222576708968190643735578161547893435503117e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.55114022296440557292153472011659421232773792851704649e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -3.29558714084043876803654097222727620781593122323144986e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.20655790191462561334686083455624901268509187198531392e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -2.29211265190736317625627743320462821391020598683096969e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.16547271274032555838193953745182811226836340887412893e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -4.04158977282328622262713955181985365200737836471868621e-06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 9.04099527601588222217706879074903022432460969180041885e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.21026124816545340223909137077288908237190834054935310e-09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.73486546143800890651761116982457576405090186184403710e-12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -2.92426051924806467732549205594664201603409972582813132e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.61340208721806452778722682494918504979435178283255061e-17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.13517996614037348212434679983945901825382670014683797e-20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 4.61930266547148210568279523485628689044657261784716591e-25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -5.78463508004708787795311937422590618909194917033730955e-31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.53568640059658459444865728826184250691297909350787934e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 6.63117406343686364701920733651409741863929155565374940e-36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -2.02296711810469513023726896410540806344003955071812851e-36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 4.05776322342042086937257880311088879198393056129493081e-37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -4.82527152907979058393765521600344474439209706209233775e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.61571352467214817635971061532522193500863154360283722e-39))
      };
      T result = 0;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (k * dz + k * k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      static const T d[34] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.27149725743552824677345467452765724462183747425698430e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.29461157973340118022513089366929246812794440622918275e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 6.12714892558448280578038020250503421879861075137953016e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.78987853703364598882843888232934449281930468883391335e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.61421402429854740596701737352164170551274801368708941e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -5.35389897687788587032725281528384185518833066402901114e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 6.02677338784735277115375432621012998338558552597357336e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -5.26946128083650646734645560602832582368091536233949979e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.62745979285730759140267927796734683888575236461558969e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.98136451866684089001192332403279162096332022571108845e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.61583350640418693864271497678443168377071451407567658e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -2.98152856924573421857895007425890467500099537012952122e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 8.18020021813066494304354702578430478815286552616366363e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.76687516020609325102781997617812558189956529947937455e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 2.97286840137911020382328323869139629696357829352821391e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -3.84037960234558373817068181823488833809395626348103188e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.73663296826451629159939759241902446349258313114071094e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -2.67102106498038196040637483025426329160024817197530409e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.35813663599770253165594840156807550072235435724240491e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -4.70970369202279764117806847443483943084782006868895857e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.05355593279987378383061490338957570599931723523376630e-06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.41032914996100466638712258076005313289749145071459453e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 1.01788232911919897659739133018267424364933824143275367e-10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -3.40766909510419913528326529229511125941552363401838743e-13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 4.21073243637136845906400015214174112719748429942427344e-16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.32283620509730411652267827642614145629073466029279697e-19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 5.38291811911040907926983774874091831690977569145002154e-24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -6.74089126429842876230950047383862474772290970041867579e-30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -1.78955023078122126762117283369038212754061385626732316e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 7.72737133764019622888917725972972138033556721960526187e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -2.35738316863436492576730648996735973692694426493195922e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 4.72855077059130427852097715321796484769785970590691282e-36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, -5.62293563001680479553446157798018614684701661703654843e-37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 168, 3.04811629504314504274190599430387363985052709041515655e-38)),
      };
      T result = 0;
      T z = dz + 2;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (z + k * z + k * k - 1);
      }
      return result;
   }

   static double g() { return 2.96640371531248092651367187500000000000000000000000000e+01; }
};

//
// Lanczos Coefficients for N=48 G=2.880805098265409469604492187500000000000000000000000000000000000e+01
// Max experimental error (with 60-digit precision arithmetic) 51eps
// Generated with compiler: Microsoft Visual C++ version 14.2 on Win32 at Oct 14 2019
// Type precision was 201 bits or 63 max_digits10
//
struct lanczos48MP : public mpl::int_<201>
{
   template <class T>
   static T lanczos_sum(const T& z)
   {
      static const T num[48] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.761757987425932419978923296640371540367427757167447418730589877e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.723233313564421930629677035555276136256253817229396631458438691e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.460052620548943146316510839385235752729444155384745952604400014e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.118620599704657143233902039524163888476114389296433891234019212e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.103553323924588863191816202847384353588419783622786374048756587e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.051624469576894078907076790635986076815810433950937821174281248e+69)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.865434054315747674202246332480484800778071304068935338977820344e+68)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.291785980379681713553231795767203835753576510251486784293089714e+68)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.073927196464385740270105346713079967925505577692095446860826790e+67)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.884317172328855613403642857232246924724496526520223674336243586e+66)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.515983058669346491005379681336434957516572863544374020968683717e+65)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.791988252541273516986153564408477102509671668999707480365384945e+64)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.645764905652320236264233988360776875326874810201273735655153182e+63)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.144135487589921315512939394666974184673239886993573956770438389e+62)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.444700846549614719681016920231266383188819427952261902403138865e+61)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.721099093953481665535866508692670759355705777392277743203856663e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.100969797434901880312682514502493221610943693861105392844971160e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.418121506159806547634040503980950792234471035467217702752406105e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.417864259432558812733518752689742288284271989351444645566759428e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.665995533734965936996397899459612023184583125575089834552055942e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.444766925649844009950058690449625999301860892596426461258095232e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.053637791492838551734963920042182131006240650838206322215619662e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.150696853422753584935226676401667305978026730065639035499393518e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.985976091763077924792684854305586783380530313659602423780141188e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.269589095786672590317833654141210781129738119237951536741077115e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.718300118825405526804849893058410300716988331091767076237827497e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.001037055130874457401651655102738871459032839441218104652569066e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.475842513986568687160423191409256650108932454810648362428602348e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.620452049086499203878684285356863241396518483154492676811559133e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.169661026157169583693125067814111812572434991018171004040405784e+40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.227918466522161929152413190031319328201533237960827483146218740e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.876388843752351291646654793076860108915313255758699513365393870e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.145947758366681136606104191450792163942386660344907590963820717e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.853323303407534484800459250019301328433169196161471441696806506e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.154628006575221227908667538321556179086649067527404327882584768e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.357526820024103486396860374714568600536209103260198100884104997e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.431529899588725297356982438015035066854198997921929156832870645e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.345565129287503320724079046959642760096964859126850291147857935e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.118851309483567225684739040233675455708538654675741148330404763e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.153371780240325463304870847387326315142505274277395976930776452e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.146212685927632120682088036018035709941745020823689824280902727e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.771109638413640784841091904266004758198074452790973613270876444e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.247775743837944205683004431867637625466576857881195465700397478e+14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.570311375510395966207715903995528566489264305503840005145629111e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.307932649387240491969419239876926639445902586258953887216911993e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.743144608535924824275750439447323876880302369055576390115394778e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.749690888961891063146468955091435916957208840312184463551812828e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.506628274631000502415765284811045253006986740609938316629929233e+00))
      };
      static const T denom[48] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 0.000000000000000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.502622159812088949850305428800254892961651752960000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.430336111272256671478593169569751383305061494947840000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.920361290698585974808779016476219830728024276336640000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.149178946896205138947217427059336370288899808821248000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.374105269656119699331051574067858017333550280343552000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.521316226597066883749849655326023294027593332332429312000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.808864152650289891915479515152146571014320216782405632000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.514810409642252571378917003183814999063638859346214912000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.583350992233550434239775839017811699814141926043903590400000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.478403249251559174520099458337662519939088809134875607040000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.848344883280695333961708798743230793633983609036568330240000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.943873277267014936040757307088314776495222166971439104000000000e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.331069721888505257142927693659482094449571844495257600000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.196124539826947758881834650235619760202156354268084224000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.723744838816127002822609734027860811982593574672547840000000000e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.205691767196054136766333529400075228162139411801728000000000000e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.517213632743192166819003098472340901249838381523200000000000000e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.571722144655713179046526371841394014407124514352640000000000000e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.359512744028577584409389641902976782871564427046400000000000000e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.949188285585060392916084953872833077002135851920000000000000000e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.452967188675463645529736303316005271151737332000000000000000000e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 9.790015208782962556675223159728484084908850744000000000000000000e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.970673071264242753610155919125826961862567840000000000000000000e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.299166890445957751586491053313346243255473500000000000000000000e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.652735578141047520337049888545244673386975000000000000000000000e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.508428802270485256066710729742536448661900000000000000000000000e+40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.093294777021479729147119238554967297499000000000000000000000000e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.155176275192359061296447275633302204250000000000000000000000000e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.907505708457079284974986712721395225000000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.195848283940498442888394846136646210000000000000000000000000000e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.305934675041764670409270520636101000000000000000000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.238840089027488915014959267151000000000000000000000000000000000e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.846167161648076059624793804150000000000000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.707826341119682695847826052600000000000000000000000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.648183019818072129964867660000000000000000000000000000000000000e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.059080011923383455919277000000000000000000000000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.490144286132397218940500000000000000000000000000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.838362455658776519186000000000000000000000000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.970532718044669378600000000000000000000000000000000000000000000e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.814183952293757550000000000000000000000000000000000000000000000e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.413370614847675000000000000000000000000000000000000000000000000e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 9.134958017031000000000000000000000000000000000000000000000000000e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.765795079100000000000000000000000000000000000000000000000000000e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.928125650000000000000000000000000000000000000000000000000000000e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.675250000000000000000000000000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.081000000000000000000000000000000000000000000000000000000000000e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.000000000000000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      static const T num[48] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.775732062655417998910881298714821053061055705608286949609421120e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.688437299644448784121592662352787426980194425446481703306505899e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.990941408817264621124181941423397180231807676408175000011574647e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 9.611362716446299768312931282360230566955098878347512701289885826e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.401071382066693821667231534775770086983519477562699643517826070e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 9.404885497858970433702192998314287586471872015950314081905843790e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.115877029354588030985670444733795075439494699793733843615128537e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.981190790128533233774351539949086864384527026303253658346042487e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.391693345003088328615594164751621620795026048184784616056424156e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.889256530644592752851605934648543064680013184446459552930302708e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.083600502252557317792851907104175947655615832167024966482957198e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.168663303100387254423547467716347840589509950430146037235024663e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.123598327107617380847613820395680616677588511868146055764672247e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 9.689997752127767317102012222013845618089045780981297513260591263e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.534390868711924145397558028431517797916157184545344400315049888e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.304302698603539256283286371502868034443493795813215278491516590e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.393109140624987047793401361048831961769792029208766436336102130e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.978018543190809154654104033779556195143800802618966016721119650e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.053360999285885098804414279382371819392475408561904784568215676e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.134477518753880004346650767299407142912151189519394755303948278e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.294423222517027804991661400849986263936601088969957809227734095e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 9.411090410120803602405769061472811786006792830932395177026805674e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.546364324011365762789375386661337991434000702963811196005801731e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.228448949533845774618310075362255075191314754073111861819975658e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.912781600174900095022672513908490962899309128877584272045832513e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.145953154225327686809754524860534768156895534588187817885425867e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.085123669861365984774838320924008647858451270384142925874188908e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.630367231261397170650842427640465271470437848007390468680241668e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.732182596346604787991836614669276692020582495778773122326853797e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.604810530255586389021528105443008249789929772232910820974558737e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.866283281281868197964883431828004811500103664332499479032936741e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.194674953754173153419535571352963617418336620849047024493757781e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.894136566262225941799684575793203365634052117390221232065529506e+22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.728530091896234109430773225830735206267902257956559214561779937e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.558479853180206010560597094150305393424259777860361999786422123e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.183799294403182487629551851184805610521945574359855930862189385e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.411871423439005125979602342436157376541872925894678545707600871e+14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.146934230284030660663814250662713645615827253848318877256260252e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.448218665084135299794121636822853382005896647323977605040284573e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.512810104228409918190743070957013357446861162954554120244345275e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.586025460907685522041021408846741988415862331430490056017676558e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.540359114012197595748944623835295064565126012703153392373623351e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.845554430040583564794301575257907183920519062724643766057340299e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.408536849955106342184570268692357634552350288861587703063273018e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.030953654039823541442226125506893371879437951634029024402619056e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.454172918244607114802676127860508419821673596398248024962237789e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.155627562127299657410444702080985966726894475302009989071093439e-09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.725246714864934496649491688787278190129598018071339049048385845e-13))
      };
      static const T denom[48] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 0.000000000000000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.502622159812088949850305428800254892961651752960000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.430336111272256671478593169569751383305061494947840000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.920361290698585974808779016476219830728024276336640000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.149178946896205138947217427059336370288899808821248000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.374105269656119699331051574067858017333550280343552000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.521316226597066883749849655326023294027593332332429312000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.808864152650289891915479515152146571014320216782405632000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.514810409642252571378917003183814999063638859346214912000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.583350992233550434239775839017811699814141926043903590400000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.478403249251559174520099458337662519939088809134875607040000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.848344883280695333961708798743230793633983609036568330240000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.943873277267014936040757307088314776495222166971439104000000000e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.331069721888505257142927693659482094449571844495257600000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.196124539826947758881834650235619760202156354268084224000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.723744838816127002822609734027860811982593574672547840000000000e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.205691767196054136766333529400075228162139411801728000000000000e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.517213632743192166819003098472340901249838381523200000000000000e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.571722144655713179046526371841394014407124514352640000000000000e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.359512744028577584409389641902976782871564427046400000000000000e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.949188285585060392916084953872833077002135851920000000000000000e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.452967188675463645529736303316005271151737332000000000000000000e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 9.790015208782962556675223159728484084908850744000000000000000000e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.970673071264242753610155919125826961862567840000000000000000000e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.299166890445957751586491053313346243255473500000000000000000000e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.652735578141047520337049888545244673386975000000000000000000000e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.508428802270485256066710729742536448661900000000000000000000000e+40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.093294777021479729147119238554967297499000000000000000000000000e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.155176275192359061296447275633302204250000000000000000000000000e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.907505708457079284974986712721395225000000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.195848283940498442888394846136646210000000000000000000000000000e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.305934675041764670409270520636101000000000000000000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.238840089027488915014959267151000000000000000000000000000000000e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.846167161648076059624793804150000000000000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.707826341119682695847826052600000000000000000000000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 6.648183019818072129964867660000000000000000000000000000000000000e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.059080011923383455919277000000000000000000000000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.490144286132397218940500000000000000000000000000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.838362455658776519186000000000000000000000000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.970532718044669378600000000000000000000000000000000000000000000e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.814183952293757550000000000000000000000000000000000000000000000e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.413370614847675000000000000000000000000000000000000000000000000e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 9.134958017031000000000000000000000000000000000000000000000000000e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.765795079100000000000000000000000000000000000000000000000000000e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.928125650000000000000000000000000000000000000000000000000000000e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.675250000000000000000000000000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.081000000000000000000000000000000000000000000000000000000000000e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.000000000000000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      static const T d[47] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.059629332377126683204423480567078764834299559082175332563440691e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.045539783916612448318159279915745234781500064405838259582295756e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.784116147862702971548198855631720823614071322755242269800139953e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.347627123899697763041970836639890836066182746484603984701614322e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.616287350264343765684251764154979472791739226517501453422663702e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -3.713882062539651653939339395399443747287004395732955159091898814e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.991169606573224259776909844091992693404451938778998047720606365e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -3.317302161605094814956529918647229867233820698992970037871348037e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.160243421312714521088457044577429625205805822189897013706603525e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.109943233027050100899811890306430189301581767622560123811853152e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.510589694767723034579229465791750718722450232983242500655372350e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.447631000703120050516586541372187152390222336990410786008441418e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.650513815713423478665128697883383003943391843803280033790640056e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -7.169833252147741984016531016457108860830636610643268300442548571e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.082891222574188256195988224106955541928146669677565424595939508e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.236107424816170540654753273736991964308279435358993150196240041e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.042295614972976540486053879488442847688158698802215145729595300e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -6.301008161384761854991230670333450694872613042265540662425668275e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.626174700692043436308812511757112824553679923076031241653340508e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -7.165638597797307942127436742547456896168876912136407736672893749e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.193760947891421842393017150194414897043594152709554867681454093e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.102566205604210639065160857917396944102487766555058309172771685e-10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.915816623470797626925445072607835810426224865943397673652473644e-13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -8.588275837841705058968991523347781566219989845111381889185487327e-16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.200550301285945062259329336559146630395284987411539369061121774e-19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -3.164333226683698411437894680594408940426530663957731548446585176e-23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.066415481671710192926882432742434212829003971627792457166443068e-28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.794259516500627365643093960688415401054083199354112116216326548e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -4.109766027021453750770079684473469373477285891593627979028234104e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 7.857040454431507009464118652247309465880198950544005451066913133e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.257636833252205356462338019252188768182918234805529456629813332e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.657386968948568677903872677704817552898314429680193647771915640e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.807368757318279512579151153998249666772948741065806312921477647e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.661046240741398691824399424582067048482718145278248186045239803e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.310274358393495831279259654715581878034928245769119610060724565e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.979289812994200254512860775692570111131240734486735844065571645e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -5.374132043246630393307108400571746261019561481928368054130159659e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.807680467889122570534300256450516518962725443297886143108832476e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.273791157694681089776609329544693948790210894828257493359951461e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.971177216154470328027539744763226999793762414262864963697237346e-36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.645869582759689501568146144102914403686604774258048281344406053e-36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.533836765077295478897031652308024155740827573708543095934776509e-37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.011071482407693628614243045457397049948479637840391111641112292e-37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.753334959707221495336088007359122169612976692723773645699626150e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -2.217604773736938924265403811396189809599754278055061061653740309e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.819104328189909539214493755590516594857915205552841395610714917e-40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -7.261124772729210946163851510369531392121538686694430629664292782e-42))
      };
      T result = 0;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (k * dz + k * k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      static const T d[47] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.201442621036266842137537764128372139686555918574926377003612763e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.185467427150643969519910927764836582205108528009141221591420898e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.424388386017623557963301151646679462091516489317860889362683594e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.527983998220780910263892115033927387104053611029099941633323011e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.966432728352315714505545454293409301356907573727621630702827634e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -4.210921746972897898551337991192707389898034825880579655985363009e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.525319492963037163576188790739239848749059077112768508582824310e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -3.761266399512640929192286468240357629226481512485264527650043412e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.449355108314973517543246836489412427594992113516547680523282212e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.258490177973741431378782429416242097479994678322390199981700552e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.114255088752286384038861754183366335220682008583459292808501983e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.641371685961506906939690430062582517060728808639566257675679493e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.139072742579462987548668350779672609568514018384674745960251434e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -8.129392978890804438983060711164783076784089453197491087525720250e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.227817717944841986447189375517242505918979312023367060292099051e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.401539292067249253713639886818857395065226008969910929456090178e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.181789081601278618540976740818676551399023595924451938057596056e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -7.144290488450459735914078985115746320918090890348935029860425141e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.977643331768050273059868974450773270172308183228656321879824795e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -8.124636941696344229278652214634921673116603924841964381194849043e-06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.353525462444406600575359080915245707387262742058104197063680358e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.250125861423094782405286690199652039727315544398975014264972834e-09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.573714720964717327652547152474097356959063887913062262865877352e-12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -9.737669879005051560419153179757554889911318336987864449783329044e-15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 4.762722217636305077074994367900679148917691897585712642440813437e-18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -3.587825185218585020252537180920386716805319681061835516115435092e-22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.209136980512837161314713015292452549173388035330975386269996826e-27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.034390508507134900778125110328032318737425888723900242108805840e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -4.659788018772143666295222723749466460348336784193790467337277007e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 8.908571128342935499766722474863105091718059244706787068658556651e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.425950044120254934054607924023969978647876123112048584684333719e-33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.879199908120536747953526966437055347446296944118172532473563579e-33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -2.049254197314637745167349860869170443784687973315125511356920644e-33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.883348910945891785870183207161008885784794173754432580579430117e-33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.485632203001929498321635338807138918181560966989477820879657556e-33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.018101439657813295290872898460623215815148336073781084176896879e-33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -6.093367832078140478972419022586567008505333455627897676553352131e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 3.183440545955440848970303491445824299419388286256245840846211512e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.444266348988579122529259208173467560400718346248315966198898381e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.636484383471871096369253024129613184534143941833907586683970329e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.866141116496477515961611479835778926021343627571438400431425496e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 5.140613382384541819628458619521408963917801187880958447868987984e-36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -1.146386132171160390143187663792496413753249459594650450672610453e-36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 1.987988898740147227778865012441676866493607979490727350027458052e-37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -2.514393298082843730831623322496784440966181704206301582735570257e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, 2.062560373914843383483799612278119836498689222815662595453851079e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 201, -8.232902310328177520527925464546117674377821202617522000849431630e-41)),
      };
      T result = 0;
      T z = dz + 2;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (z + k * z + k * k - 1);
      }
      return result;
   }

   static double g() { return 2.880805098265409469604492187500000000000000000000000000000000000e+01; }
};
//
// Lanczos Coefficients for N=49 G=3.2804746093749997726263245567679405212402343750000000000000000000000000000e+01
// Max experimental error (with 60-digit precision arithmetic) 88eps
// Generated with compiler: Microsoft Visual C++ version 14.2 on Win32 at Oct 14 2019
// Type precision was 234 bits or 73 max_digits10
//
struct lanczos49MP : public mpl::int_<234>
{
   template <class T>
   static T lanczos_sum(const T& z)
   {
      static const T num[49] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.5742499367008530642617319352897520430205552639474773193534726699460038382e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.2055424563420473601452029826557352738661904429902299059433513914586136617e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.5124573975150897396579361228988665521507891014129938036244581220575769357e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.7657855450073472861630911067774749012357295936579551712546801770340228765e+73)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.2200724511833181761435003548793565030399736349888521665887428093215263908e+73)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.6969626159645181170657767707096364240775930310958541253756412108457464465e+72)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1902838064351166646672468793217979804852505454197638767497901741630038659e+72)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.0815524946284433112181691476757266928589723434617417902816279882371678020e+71)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.1085775474754412572785597234367295850853434122501232697548485520290211685e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.0248876219420480494671692552535249842631934230477107608175338135254767453e+69)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.5717898581583284547576613790182065259217990496576104453319232531339965508e+68)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.5987266955225124754687270016251985476100549398788271466921127242855554305e+67)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.1277201593916609119851488898621364059725185338025069187337522019022532875e+66)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.3266912844842724918464413240523721967983751112927460422203771012324285814e+65)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.4198349306520385923755583127019190414798256355752152907538684107889754415e+64)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.5954968135971263467860541616149005468520952700568675174669869313790120813e+63)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.5697573028906337479684490638442403255222113759866980140004921288934014321e+61)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.2372180596212313667398562038281895205587396825353276904411403952354592724e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.6216537416683580103634264243351068913266542152026155538733152213171424662e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2028581946420676660634538276758371308595650705510961727364638531545389749e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.0668609223671045472881034712995505846341571161239968651357297157517601241e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.9620881267196496193861934343921326919697524102893094861578977849944522047e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.9917155089797622879653341069010820156774380191279115491781432309218621096e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.2942244374439765758128985411022615606270110470067503605960457871147797553e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.9351028342882568351470704633625779062510316389795468890729515726996746376e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.9314927110707228066333377498237626208657080992885612416826068446390297282e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.9556005176089621449765653988254594966772984328873034691437246673533838584e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1708002669765411191695114973458120651885979389757503801447241827647515479e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.5453624823926493383313636717841737625222638948016893410076575426986147237e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.0870291142516923896770486348199039577071847823857246263136765980806245584e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.3337301613565742543387644582305705841930917569470575861418313137426853827e+40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.5696392889303372910143506267477452989689410812923616782357472172107995484e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.4143864855807616500285749805611035801672647319221651520900664580210701417e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.3883950248556863481116504212968387826172322376709179509012607292972789664e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.3257340722619328574784440823494550029955992275874284892851081786012219555e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.0054658629154523660536922265411423770787857332470130847913952020455172951e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.2273347597734780854875327656763513602785317406227796368685036961474823856e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.9014441058506987415445336961216002315976811002826148787399813116200761040e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.1007719003538290551470147047960192204273669359908829249030378418530058831e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.0380905019192393560853564871253022486270254732200016399586810852416299561e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.9744635849298784873006356360567530952404226257040467086536310864153303851e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1124649237555162417911032704785227088470234269296071198905173077670702197e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.3521744565441850034376761928380748678988834367585778746740882700228177944e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.1551192264984335079683102449711150606277601324509245213720618349409077054e+14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 7.0649463728846790219460420909345995584550543072052058591106589782466481621e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.8110794685236845669275412185521099716923991533427218655037515888279356458e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.4051945147325155571782612757500183884572379142638708146516926708747795902e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.1761279721833493440329993975995822357104129877086617643903445974526391316e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.5066282746310005024157652848110452530069867406099383166299235756022800280e+00))
      };
      static const T denom[49] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 0.0000000000000000000000000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.5862324151116818064296435515361197996919763238912000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1477605944577727245447890951265834050463405543784448000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.3368731677410579748749120694395208342752220248276992000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.9393177179482022750532799808826502923430631529093529600000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.5873212662073383100750664140824866318496576298496819200000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.7087596791971826323557398537439095283663043689996772966400000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.8537931401160691803777386867476912131700643521105494016000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.7128473077968876977396389430116077066613422855709615718400000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2893230704461912298064838143702096489032830938340968366080000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.7731846263715878554484243293204825543527859328977818943488000000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.4350612763444239870720412999269509820736318433853587128320000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.0384549286435665533353268142058310243161527793802332119040000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.8399900970142989644612517467287880620408209836099149824000000000000000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.3548923093755049924589156254733610823950920495095216128000000000000000000e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.2977252822627446721481004001665655765203461552290590720000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2090496144637581445624377322208214384344648810140669440000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.4036595841089057320815648092220077464036379804960768000000000000000000000e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.4604307712625044108337677046126892768963323598980608000000000000000000000e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.3661432041590027825770657688785384893903477321470720000000000000000000000e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1520697686278361431114988925105292244781602931070400000000000000000000000e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.7781340723597395269058455794580578514153013123200000000000000000000000000e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.0542743368034560471670911883883927910588971816800000000000000000000000000e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.7852178643724903498642955979619870805662919592000000000000000000000000000e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.1476757456360244186066663869698554305163293290000000000000000000000000000e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1067024107708881097170625529475996208174256000000000000000000000000000000e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.1816971152081755906884039315242368042580680000000000000000000000000000000e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.2046914254271439983058171150950882746907200000000000000000000000000000000e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.5226232704255673172404214340314876574740000000000000000000000000000000000e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.9917039581671863252346910306123579600000000000000000000000000000000000000e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.5279926429774219665504424895636324120000000000000000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.7496375812101278379807519908356136800000000000000000000000000000000000000e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 7.1781895168846844604663013761970710000000000000000000000000000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.6915825748773446395251490146656000000000000000000000000000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.5888455419743269266732720488720000000000000000000000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.8324723604341765969313138528000000000000000000000000000000000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1625859075857974372785469560000000000000000000000000000000000000000000000e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.7594478264056101488213120000000000000000000000000000000000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.3541746402920221829579200000000000000000000000000000000000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.7645128331397711271280000000000000000000000000000000000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.8231991756227354271000000000000000000000000000000000000000000000000000000e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.4784681412721648000000000000000000000000000000000000000000000000000000000e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.8427136416481320000000000000000000000000000000000000000000000000000000000e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1374881704208000000000000000000000000000000000000000000000000000000000000e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.6720141346000000000000000000000000000000000000000000000000000000000000000e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.1948624000000000000000000000000000000000000000000000000000000000000000000e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.1833200000000000000000000000000000000000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1280000000000000000000000000000000000000000000000000000000000000000000000e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.0000000000000000000000000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      static const T num[49] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.9156595384723483714761513394165493705512008185104458104037189919966679583e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2490942626049052041642905803851633038401948356209065994820093902482046488e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.5657016133969059292054137486140786946523144967749946090276991389460410256e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.8317575261281828541202228096164333127703313535243593976139496326255603885e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2573232282907489549580534670684496645767812633106925447384772399512629009e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.2264367876545276188760221330510623567338280567044843125670581811561508361e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.7410929642927466039608438519399985289266866350365951941371580436768272283e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1788733746089745519939950216688781077254011132768568000829003039800512115e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.7605221646260701306188822826505483714178735327624200723152581930359377597e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.2794682649344687478785548483237660272987192396957769664741689026760236827e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.5892026994266264487972781655291611074377224857532819461836222961458821278e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.6044581976409422084009086165534449955654966215396670231211076839220955529e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.3377067867856641370165698166594859838986332949917530009978864573093452852e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.8840494251009843798733373985085724953611000664053469151896413377742842585e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.3704573764322216880093728898018690177362850852138123977146826586562660042e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.0359898089375286965109471568127438755772335943188122764617864358438172366e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.4197682330664983610360781372995395147813463692379312665516728111546651975e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.9660635239520142135892297007271080519825235074724441273950199516718744758e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.4847561142331378363753305417832438611221296838300212519621989387291450611e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.8123071733863185211372790750342775813125583218947155817719674827525929182e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.8695828952857913310728321803468643286198192850402657393713980287261259638e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1112155264857128661506293440611411838336453876178977662636328219208269490e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.9597114546218224309028736968589856884375435031334606638164700453351137400e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2993187112308342042984815542969104496107563160266146124117192158716725839e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.9276492438289459104440745719026053099029450534209231660211458499968839130e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.0938880168568707970634270396301729019623679984385010064117584299386152912e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.8065713069852021979452626742757704686014267992443726831364663812595589478e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.6307492378187334600766084298235996134621924808224715558467640454423842655e+31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.4415490682866089784509877074414842364885644695603000622810515902853301434e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.8810050162690958474575145223002401297790616292489104417882204295223485358e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.2860958353933214242194744963395870718749107497800918932121607454081905640e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.8895474422826652624342608217034044478811180021465482602156553101973142851e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.3673716858987747316948988327383411119224741422002091782489841976044254452e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.9189949269921720477131833122497011042624119182630703457459218878350720360e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.4498506458943304056122167501858896499220570899589497703614276271780187559e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.8348122127750881334333964077955476539118633228361300065573128691903485670e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.9604661829895068017545272035886457795302072787311442846001622570781400887e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.7759001843251065234667097315555616424286286440771094985162292676933917344e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.3224448199826427582356965293979283708473565146649764333289536983716495991e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.7206022964144885856033111994076522875074747824871090843993137657458226608e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1182242847179807470291498763391918222721736346452719834560856782062915656e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.3003709115482840946768994349442745995208542698397935766566883071074706882e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.0311683127686139178455785110448360190741004055855567244689436514834230540e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2205374026276431316033131273698903101132908812758212963153267146702515001e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.0011852660581895977871209152278420195617235126960126683426738428884540341e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.0256927799097601645189517227865464454408483715314872385869625463116763013e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.9285092060573954544191464750091257510108951419853791428165851425157303961e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.3651222287552713546296738256776548189076669986821606599921179638781057448e-11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.4196122080178359418232367999588103764109880361922878363699332724809241837e-14))
      };
      static const T denom[49] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 0.0000000000000000000000000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.5862324151116818064296435515361197996919763238912000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1477605944577727245447890951265834050463405543784448000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.3368731677410579748749120694395208342752220248276992000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.9393177179482022750532799808826502923430631529093529600000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.5873212662073383100750664140824866318496576298496819200000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.7087596791971826323557398537439095283663043689996772966400000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.8537931401160691803777386867476912131700643521105494016000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.7128473077968876977396389430116077066613422855709615718400000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2893230704461912298064838143702096489032830938340968366080000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.7731846263715878554484243293204825543527859328977818943488000000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.4350612763444239870720412999269509820736318433853587128320000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.0384549286435665533353268142058310243161527793802332119040000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.8399900970142989644612517467287880620408209836099149824000000000000000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.3548923093755049924589156254733610823950920495095216128000000000000000000e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.2977252822627446721481004001665655765203461552290590720000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2090496144637581445624377322208214384344648810140669440000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.4036595841089057320815648092220077464036379804960768000000000000000000000e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.4604307712625044108337677046126892768963323598980608000000000000000000000e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.3661432041590027825770657688785384893903477321470720000000000000000000000e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1520697686278361431114988925105292244781602931070400000000000000000000000e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.7781340723597395269058455794580578514153013123200000000000000000000000000e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.0542743368034560471670911883883927910588971816800000000000000000000000000e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.7852178643724903498642955979619870805662919592000000000000000000000000000e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.1476757456360244186066663869698554305163293290000000000000000000000000000e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1067024107708881097170625529475996208174256000000000000000000000000000000e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.1816971152081755906884039315242368042580680000000000000000000000000000000e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.2046914254271439983058171150950882746907200000000000000000000000000000000e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.5226232704255673172404214340314876574740000000000000000000000000000000000e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.9917039581671863252346910306123579600000000000000000000000000000000000000e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.5279926429774219665504424895636324120000000000000000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.7496375812101278379807519908356136800000000000000000000000000000000000000e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 7.1781895168846844604663013761970710000000000000000000000000000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.6915825748773446395251490146656000000000000000000000000000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.5888455419743269266732720488720000000000000000000000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.8324723604341765969313138528000000000000000000000000000000000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1625859075857974372785469560000000000000000000000000000000000000000000000e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.7594478264056101488213120000000000000000000000000000000000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.3541746402920221829579200000000000000000000000000000000000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.7645128331397711271280000000000000000000000000000000000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.8231991756227354271000000000000000000000000000000000000000000000000000000e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.4784681412721648000000000000000000000000000000000000000000000000000000000e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.8427136416481320000000000000000000000000000000000000000000000000000000000e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1374881704208000000000000000000000000000000000000000000000000000000000000e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.6720141346000000000000000000000000000000000000000000000000000000000000000e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.1948624000000000000000000000000000000000000000000000000000000000000000000e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.1833200000000000000000000000000000000000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1280000000000000000000000000000000000000000000000000000000000000000000000e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.0000000000000000000000000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      static const T d[48] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2066789802575522178795921101394468758895980798225454578416417743837821182e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.3679788352265125219923749199991620615257401271523573821302279091870986385e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 7.2630033804072373059579276357751429378169881873880556882995810292143631169e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -2.3999752086438541157776701982236284193725611536088710297195418612889088941e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.5326869500703014940743489551921365917396663139743104045090175695793718397e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -9.4539806536100390983482685579769998888148275478188901957443850396223701663e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2418837249375022042463276831798225069762119451024109889652570545798197307e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.2837532371348753406818677898802253115337338144774680681639085300947433821e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.0603676137809014970843333435453530863140601935101507289758344098164636049e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -7.0675156715558786066521224677848991893948220792765278914746791001795702147e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.8231576762102015440131857562883207655276091648820750106248419664382342083e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.6828357670521487851037941818117892204609675060095269337474530867886953446e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 6.0263919393474265394410446439494942658409182482953883911882462037294240391e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.7511174204844870012138672336839720481677091590654128683087419414952818657e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.1077240937707245791059934044373222445269304780721564415689327491376493334e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -7.7197427650993388314150434222790763404096237849031957731092258677118877326e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1502942251829477374710555866529647619927300737021040120015632037419548655e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.3407293979752550880132165447474755706150471817851620226242890535888257338e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2014334746887994911502478644617408322636134206225589060774099980181312765e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -8.0983646402105726518672426329015155996093741444426226630025345020847577111e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.9942365899792195557178345302806016444352718230344395806792868608352214871e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.3916756170744060855799456788870481671060079696580348774412746337324507629e-06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.2738029264164285805306112157985005833945884241358728352406636162060451936e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -4.9005795832591317070122324350352498562565929180970558213118848371143071658e-10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.3120237284050523317140903263269508807725654665206264201337566990649593815e-12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.9993164488488861127293338637233872747923249255635013388317575640760246175e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.1770960183589659510713725797228496911736699917131125983799031646762143273e-17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -3.1067143632180728240506803448567684891294208768806045480286272494772691950e-20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.6244126131167038529019905678720925754781858101865433373489388749513540167e-24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.2602412034992451089948288841723108481473094784377158526788906759478344724e-28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 8.5993275050503412644261961442095022244184843038273154020407386192200117900e-35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2376677581089123203915935849311668028430409001498177367770610457728550321e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.8150848887085059253468109780286976776338724012461857375172227849044365287e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.2888861291942792803542790121202831584212092799440963290014968720997192671e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -2.4367686092860277330824120317504561350242498766190734227637416697402218133e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.1658834966582889258321173900424691594425653659844978485212529673045052341e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.6078489664179099912530462853749930668959707942650569176553871936075352345e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.0044570187922899426305396174281740469543589839974788550682447627990764739e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -5.3214978442929458907632954943316397146183013077035413667972521651417588659e-40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.3990788346383540827394552555896039808285381860296208830197355311239959128e-40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -9.1813962257175782149431101470253920074573257556345440519764098954533396958e-41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.9566302522237219982311281591833653462640003193980940062047608523626590962e-41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -7.8844939069206979681852040267695059179818896523029249336074241598129744909e-42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.6980186103098555864348162964443149317288303456071195638226831779183346735e-42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -2.8409438292612321309493446423416582883045091540877596567806910261300252610e-43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.4666144269083503254700296220598377304876489203081567665937809718753144311e-44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -2.7449599736044707057960973836294763832382889203992176649603299199814964095e-45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.0585685375038567471031307628591396658405283305367035292669018178549132199e-46))
      };
      T result = 0;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (k * dz + k * k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      static const T d[48] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.5455231550124674608383465879384903101217726809855880053735498558660305251e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.7521171744934997797979463716947798398649024826461953154035391579443199892e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 9.3025071978608000907597512615615521968983853964414599121560447447486249179e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -3.0739055847506874449724863962235222139048946552417933333489374601444889442e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 7.0863054140081769874319573695951392605490718603916106954837774653056576623e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.2108726717088904191265369690298547113580473039391314603986085479556111007e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.5906136463190712486212821488571957306185167770225864012898706881699647057e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.6442404200085453327213398272100872399623441451383083438265316066112956103e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.3581264998698427315063870744999865373965234956609066658837559718999590062e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -9.0521251281527290170664891261612911437995934455165647349853355253750334406e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.8967279703383588010975605180679227794424995039989024281878608594596616855e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -2.1553882073152030175586351692198830408611546985761579133373632331779244638e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 7.7186463308194716629084448574149234823489365098510570780593191406061889318e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -2.2428438423008152181399734545002360362921592139304278927551479335858821905e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.2612026936693838529789647834479279275034238291000514930718615542291190074e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -9.8875023012785484029602598574637362767234771045911074782572204578687659312e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.4733051533871237880957666174363507890599952187881138769651021452363258297e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.7172158984110335489611392604457344891052843492969255546784144895328933329e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.5388046735862609465152559757953131111281278387717331973511139237319990770e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.0372443934101018493930850177650436102881381002613480220638309361418767946e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.1158470789748281012838916430380981653860595791320888251473386188428788856e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.7824681838708087799263224394519010740357745634535755030787340328399772549e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.1931104382411821412587448903202613576618187938526144451753089437172395491e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -6.2766977322268088075311831528842893348256386710249195361734612266778385323e-09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.5228711415779960574453455948991424203497366984285083000319325448210101974e-11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -2.5607389508299290889393825247382548268608977769597750095981719736417860627e-13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 5.3500547558276404981577111880189245551703390604657164900460067017086475257e-16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -3.9791021994420737378470356651659599587167405408055854179971153565586792084e-19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 7.2037883058679724499875763987161800762132311085398211645162684109748694340e-23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.6141260374761318937819974670878746263129124551434880143657761001893177100e-27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.1014080790364138694339981796904397090336135156681914457153169895058123962e-33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.5852138055486955162695207029733124482316515915389725207819834510328826882e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -2.3247738377057606447050864939101938086522392428747461590714144146282468974e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.9316218892795919921760137631186915684395320376859960566984260379803979907e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -3.1210308380902232855988541317886447706766253774127605015790322866153349557e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.7740792289514179382145837836403663560484901236914958039358936821942590582e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -2.0593445713551373294655454085927889987266290706104795764516423812212247302e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.2865158059079907658671256163529513022184187924172645112729579745062723961e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -6.8158128816897517506937531950420472376185024126168447501964316612370015014e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.0727575024490003112672747176557026577347452778581212406977128868690751231e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.1759598612683119428775753540097890056929269067220807101660487076215190882e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 3.7868733858667171123449703661037140974108364228765817302050018596960416788e-40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -1.0098516753892371349956889909483048597911093401064123370095186647248754503e-40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 2.1748345026411581847490576071466694584307900894689508331424519370338054748e-41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -3.6387013796128814207293408288388785778364579569320260571422484577307186443e-42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 4.4400648009494470795051094944915668879009266928618607576647564429190022802e-43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, -3.5157645638963794722990984515235056341682398223940346512667816730324057574e-44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 234, 1.3558222299776027787073219176369904651832196778690432544150792644775402371e-45)),
      };
      T result = 0;
      T z = dz + 2;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (z + k * z + k * k - 1);
      }
      return result;
   }

   static double g() { return 3.2804746093749997726263245567679405212402343750000000000000000000000000000e+01; }
};
//
// Lanczos Coefficients for N=49 G=3.71521093750000019895196601282805204391479492187500000000000000000000000000000000000e+01
// Max experimental error (with 80-digit precision arithmetic) 75eps
// Generated with compiler: Microsoft Visual C++ version 14.2 on Win32 at Oct 14 2019
// Type precision was 267 bits or 83 max_digits10
//
struct lanczos49MP_2 : public mpl::int_<267>
{
   template <class T>
   static T lanczos_sum(const T& z)
   {
      static const T num[49] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.29576393886533812647621748658593518049105182771886345672193119325886860615058992916e+76)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.64940250521075704286553287922703978233859608760738220959805633167660346858918484441e+76)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.02782279850582131264593341856583603623569135771086732311690891701018783486399941058e+76)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.17871318833644489595452604622172623859528919270761702889241650479860350523603645576e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.24636895901125366174317752785594004344350462222105593541765135970641727807408046603e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.90766536403132277814150094353285923029162708480029761430241108404154368124040257524e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.52382781721986299602475060387708982171655203328511435511051797995280492002324724681e+73)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.78482374388380268613213967834073614067861524793728730151201542394918426634868617345e+72)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.19325027080949586393337406155391099309209314019751861163367432285571050834050337943e+72)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.40545316209557907826875120370637029196186181048557227564363347693674382156389673345e+71)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.45248231502915687322615277256109191586795311001308546524559632502753988355826682144e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.32951710710760475384871665492267398285929602882057694753684715739127717537388087360e+69)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.08609708232067139682194681883946299264352243218047148918541516101215386723695461807e+68)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.96790018517492120048938623906554309076682107791064003018940091271134763465693455777e+66)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.27669033307366420506859833929956397578574575080825866111910129757118415016658467329e+65)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.16802337748974823514819525636942825363065912689755393017268041605495575074123213602e+64)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.73054387508326490801862749035024819220968371854706600225017512476962193917082134622e+63)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.62670382432299595501509896985706425528571034019439577742070022359804568860138431456e+61)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.93419603580003232059946064747952976400826625336761312144282727389088893514776127326e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.64477326525726493380322761443767318183513031346169727378679051953394311858102484110e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.31417565379195208375864005040359979626430954027274059125918083665481151874801228147e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.22873281009357235993917586746434823288058696112725892939606733669981526181084161579e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.24038151344262567906410045463097784963500265667230009248051943941184141216283364012e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.16635103397331843098281133679963613292528445793288854356269441338108669168034969473e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.97226028792130808404381304085761534201244857265870272466382244177818968609788744866e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.51722539705401822010051258989781033530766592277336923172687967760117181597042914505e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.55143196851124103497589961680541188447723645252871343011474640057780989215029733411e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.65631576916353367608977102104279419518759334245125488959476099210969778356845251386e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.51913643147679636800192101181821796326331285929097206716389984328795132228290343736e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.77142070954017798440735921701401232362942654167866941640116269309928860707230400223e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.64265304588830405590756890757870027959619641212378439580733250347485559404826031214e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.12961287129418300922993605074391681508051626956297216905205033442654279681009631622e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.00163852645278376379452857266851549886004284453902427442045925921205957134756527439e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.28416298426892934885339736163488241290345416925540130172490965458139411097784230624e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.49793255972439978183906387866659603599388803095306163216192230005161964280456557158e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.58404837985235629013707612365235872533571223333496498253968477010731743022987544752e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.51210388163173191310089232642285522934883375636207281562370177051035472857070373332e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.29625175202138639315345810129285907904834629385738665013524881366645746417930602439e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.91704411637684940906309546543695853685028152118727394684725500498007979899474115290e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.71980627815396304296720837491060731114276254036666247979635295719640774180255608217e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.99515302992882524356421852289147703799628458429164679118485650618041455793969699036e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.05962638033903929635201902285319054754141795004681756534441119215276710860621598251e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.06860641169958199340995275230206293988239276251216776246763641876943705626715807705e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.34257282855413383182163862298720605042203946165367698145431759335580579509759239742e+14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.00325103633255410612780052993328410422609156975902055798191507344081921284741604327e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.35516908008075882088679988477934347017452660082246021362543874768562840957846903854e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.05605466237751558310956675128455566243720955874966620350919363272463056561354608608e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.55729572460925281137543520027396380407277254872226907448088245735142585109299314309e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.50662827463100050241576528481104525300698674060993831662992357634350984153163470001e+00))
      };
      static const T denom[49] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.58623241511168180642964355153611979969197632389120000000000000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.14776059445777272454478909512658340504634055437844480000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.33687316774105797487491206943952083427522202482769920000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.93931771794820227505327998088265029234306315290935296000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.58732126620733831007506641408248663184965762984968192000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.70875967919718263235573985374390952836630436899967729664000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.85379314011606918037773868674769121317006435211054940160000000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.71284730779688769773963894301160770666134228557096157184000000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.28932307044619122980648381437020964890328309383409683660800000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.77318462637158785544842432932048255435278593289778189434880000000000000000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.43506127634442398707204129992695098207363184338535871283200000000000000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.03845492864356655333532681420583102431615277938023321190400000000000000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.83999009701429896446125174672878806204082098360991498240000000000000000000000000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.35489230937550499245891562547336108239509204950952161280000000000000000000000000000e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.29772528226274467214810040016656557652034615522905907200000000000000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.20904961446375814456243773222082143843446488101406694400000000000000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.40365958410890573208156480922200774640363798049607680000000000000000000000000000000e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.46043077126250441083376770461268927689633235989806080000000000000000000000000000000e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.36614320415900278257706576887853848939034773214707200000000000000000000000000000000e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.15206976862783614311149889251052922447816029310704000000000000000000000000000000000e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.77813407235973952690584557945805785141530131232000000000000000000000000000000000000e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.05427433680345604716709118838839279105889718168000000000000000000000000000000000000e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.78521786437249034986429559796198708056629195920000000000000000000000000000000000000e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.14767574563602441860666638696985543051632932900000000000000000000000000000000000000e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.10670241077088810971706255294759962081742560000000000000000000000000000000000000000e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.18169711520817559068840393152423680425806800000000000000000000000000000000000000000e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.20469142542714399830581711509508827469072000000000000000000000000000000000000000000e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.52262327042556731724042143403148765747400000000000000000000000000000000000000000000e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.99170395816718632523469103061235796000000000000000000000000000000000000000000000000e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.52799264297742196655044248956363241200000000000000000000000000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.74963758121012783798075199083561368000000000000000000000000000000000000000000000000e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.17818951688468446046630137619707100000000000000000000000000000000000000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.69158257487734463952514901466560000000000000000000000000000000000000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.58884554197432692667327204887200000000000000000000000000000000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.83247236043417659693131385280000000000000000000000000000000000000000000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.16258590758579743727854695600000000000000000000000000000000000000000000000000000000e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.75944782640561014882131200000000000000000000000000000000000000000000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.35417464029202218295792000000000000000000000000000000000000000000000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.76451283313977112712800000000000000000000000000000000000000000000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.82319917562273542710000000000000000000000000000000000000000000000000000000000000000e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.47846814127216480000000000000000000000000000000000000000000000000000000000000000000e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.84271364164813200000000000000000000000000000000000000000000000000000000000000000000e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.13748817042080000000000000000000000000000000000000000000000000000000000000000000000e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.67201413460000000000000000000000000000000000000000000000000000000000000000000000000e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.19486240000000000000000000000000000000000000000000000000000000000000000000000000000e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.18332000000000000000000000000000000000000000000000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.12800000000000000000000000000000000000000000000000000000000000000000000000000000000e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      static const T num[49] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.49663610130791477461794140374370924490696733744664264226063545477852907120064438588e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.20884482942845999689198014205795770851399772847436497410154319183860398480228326286e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.53289916571146966947513241858798999095967787802671792160665282347073321227016292403e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.06257315326410705497495763375336857281692495321204603919284336971479478405381409540e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.13462097275260056500457419367381935327442299676685171001942731892832408882413435566e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.13102394952914002112762425649207244062991186339070382885956481064589709031894816496e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.04840581629048620249901745346405131124091439055732548308528857701167494671648257347e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.43838524961944675350746734141447316445830054880306039896637689653198992735066081474e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.74531483688910159100344151068655208511626316424159676205726034341092475816888011372e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.03005469109920428385867124264753643595038929217266678155910599078268836953793159656e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.06452229265586101928934888140259774581544836837635812825070642317924068191237940651e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.74401260751298616807064335204160985098707747263556384563576228484205922370857356950e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.95999059097413995107763881990067736638605879036860399501004480843044798852321320153e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.83966309607370356019329894907872727851559431468853382195539375506355518522739451561e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.86727909378079977412627692901220876068535926169729765468023963613209750668440883454e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.32183997980385285815563935626564558737112327901896765866420315909089479889407291483e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.26831322790199691523763382147990685271332413931335425536077401693909941536685446004e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.32249937786475248948857177694117967417864070990756190329454700663693499608768713815e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.88336686819032307424839603734550094343459546588744595045643380205822126447806151375e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.20545206582813581050133859680698140274276251798965738445951215893805514854794850980e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.62765065960304045423639870153223988753806599339881709607153712409175486794708548715e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.63343521375027791107490012004281169470903074417529835343633992756599689706467597925e+40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.30646566133119269769508008116142168214781437395136629206311256865360022042754477332e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.58771569575797497748917749020901655508364490119697350339089287796185236358536233400e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.37706135782308791418825498038425839394113136784263252755599976807797063043875082698e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.11197240850069894333209760920904406397403024207197696391699562264322296386093323413e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.60283960927607800898139439849090533043340726395143638152077844277751805149166560973e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.61130330576422804633963440554474730419413739922460725405112412454286488484955915987e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.11337300300297295349076762512699923957438723580465440290334845396273762855836733517e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.03117042948259288566605349078528499682057378893001882806984505986948136392761267006e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.40259403730882930955675549420130062238621449273398874217078083754343032954759973357e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.22528347572098799721775967384397848588784390912022136506439711506219967445623919686e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.34099499566403754037791880900467520578107411723617564672664393961433945628291065109e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.41161286449337414692992004461766142962293572495919749372515938105183883550234862604e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.09783271453440804906045665973485276329744852249747400433246982560607811629926185500e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.16094688076417849092104554198103013111536406662498981965117749909832209236913802539e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.10821885688578986410930893219376699456504867926107715594671444282949767588694287037e+14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.50021127722497072888981489155684575207215375480643952454970468046331143014884622646e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.26818800470068614566377728536610325478326870447541587095650095821101649087210816040e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.92493678677258766132448688907920434402513619957710479626424529957236460471189038381e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.92804216541904692956935385209461325937556277824668580633404375808334363503194457710e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.50949734377247428538055873924894146321959030826140132714559567319426159642629912422e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.64636917688205569668044105440719648394808067896062765005278597960400006426062506946e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.44976703261974330359141562652121211268881081847095890846270752149427888286399621821e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.35281305841345486358892992557075736755186129730528304574018484234880523171404766251e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.72610018227274377132777187985570646349538821991325846732587180123971090415405994340e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.97267688814850268231161594533176277702963907108885229111755006417082297177116403628e-10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.34003577384307789433517516041141801840291757339440089597599691127216013708701908122e-13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.83710441782049075320423070384744844199992659460213011403809803120117400852034616312e-16))
      };
      static const T denom[49] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.58623241511168180642964355153611979969197632389120000000000000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.14776059445777272454478909512658340504634055437844480000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.33687316774105797487491206943952083427522202482769920000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.93931771794820227505327998088265029234306315290935296000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.58732126620733831007506641408248663184965762984968192000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.70875967919718263235573985374390952836630436899967729664000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.85379314011606918037773868674769121317006435211054940160000000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.71284730779688769773963894301160770666134228557096157184000000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.28932307044619122980648381437020964890328309383409683660800000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.77318462637158785544842432932048255435278593289778189434880000000000000000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.43506127634442398707204129992695098207363184338535871283200000000000000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.03845492864356655333532681420583102431615277938023321190400000000000000000000000000e+57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.83999009701429896446125174672878806204082098360991498240000000000000000000000000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.35489230937550499245891562547336108239509204950952161280000000000000000000000000000e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.29772528226274467214810040016656557652034615522905907200000000000000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.20904961446375814456243773222082143843446488101406694400000000000000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.40365958410890573208156480922200774640363798049607680000000000000000000000000000000e+53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.46043077126250441083376770461268927689633235989806080000000000000000000000000000000e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.36614320415900278257706576887853848939034773214707200000000000000000000000000000000e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.15206976862783614311149889251052922447816029310704000000000000000000000000000000000e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.77813407235973952690584557945805785141530131232000000000000000000000000000000000000e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.05427433680345604716709118838839279105889718168000000000000000000000000000000000000e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.78521786437249034986429559796198708056629195920000000000000000000000000000000000000e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.14767574563602441860666638696985543051632932900000000000000000000000000000000000000e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.10670241077088810971706255294759962081742560000000000000000000000000000000000000000e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.18169711520817559068840393152423680425806800000000000000000000000000000000000000000e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.20469142542714399830581711509508827469072000000000000000000000000000000000000000000e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.52262327042556731724042143403148765747400000000000000000000000000000000000000000000e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.99170395816718632523469103061235796000000000000000000000000000000000000000000000000e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.52799264297742196655044248956363241200000000000000000000000000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.74963758121012783798075199083561368000000000000000000000000000000000000000000000000e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.17818951688468446046630137619707100000000000000000000000000000000000000000000000000e+33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.69158257487734463952514901466560000000000000000000000000000000000000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.58884554197432692667327204887200000000000000000000000000000000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.83247236043417659693131385280000000000000000000000000000000000000000000000000000000e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.16258590758579743727854695600000000000000000000000000000000000000000000000000000000e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.75944782640561014882131200000000000000000000000000000000000000000000000000000000000e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.35417464029202218295792000000000000000000000000000000000000000000000000000000000000e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.76451283313977112712800000000000000000000000000000000000000000000000000000000000000e+21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.82319917562273542710000000000000000000000000000000000000000000000000000000000000000e+19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.47846814127216480000000000000000000000000000000000000000000000000000000000000000000e+17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.84271364164813200000000000000000000000000000000000000000000000000000000000000000000e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.13748817042080000000000000000000000000000000000000000000000000000000000000000000000e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.67201413460000000000000000000000000000000000000000000000000000000000000000000000000e+10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.19486240000000000000000000000000000000000000000000000000000000000000000000000000000e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.18332000000000000000000000000000000000000000000000000000000000000000000000000000000e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.12800000000000000000000000000000000000000000000000000000000000000000000000000000000e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      static const T d[48] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.36662594316736437928923771651218082870790219434510733129873815200968646566980996339e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.76780165965526893166944674250172218113345269194305352999222919039854748751783233101e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.07970572889193090835848021094044050788957851893559287421136980839960026782217427894e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -4.14103116540076985021820709941614153683990156820041738016476143763552291310864305462e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.11895235901890827822967425840835669339942717903640916537446753898873381005535850448e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -2.26548471357135814759065586133601745345092666805731605105989839308832305900191949768e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.56856136042152636864745499626803467747701570337645969278665918810296351408741364840e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -4.48255074360054688238131834296168065172557135920905358164352154448630969015187687499e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.56619966800988637898793129440851486107980873809024863561425305126874295268067178801e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -3.81597333294282872581583559942700103331926900424172914177106074836091709159865752408e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.63678053487724707751260697626783885891593791493016491379265204946777211045670615960e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.51398164257933964710898171665986511410737733219980979417140085562924634596875563017e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.24278345602386851676210506047826420822246500559376040889647916551285989221630819711e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -2.88891626830459350723218099328412217578914954601038686119863249895652809029737934745e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.59719898181036663575857261276179540194083274881058286618169370083758218641306764946e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -2.64811367543293387988966817387441229757530457146808054159891759802567729960811992548e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.04201418909362106976467138015064536194810537117201428954061434907295912679183849448e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.13288541649508118041751198270788445621639048290878149127604057823054487014510567347e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.73151244996618239257199112314250678525910059470480205369884034575880974697338791461e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -2.13522814977566175652818016173551223075398445721839956765100980594656777696996995600e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.09751773607848160413869110673825438505378131216432708718954295924780177394997267064e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.61580494932170317026139049544442892885812277792439584113117496484172419053879898708e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 9.57433012294978946531397728731840486773642836205104228434646630602626154207102575497e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -4.26115172606869205202745132493304748075964471479903642930011237303391162925986198072e-06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.38298958800694421458738651547755286166537642575090914169204625632511624488640442440e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -3.15434025288605144031976010596974001946974995609543721535748916319261879201397044643e-09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.82356872364336815146383574702913568837734092081047520230601339834144374864204220536e-11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -4.65347968332361614985780283973809979932690029695650403412722596287871272868944164634e-13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.61330972151984338818972320382872350752985962075827056997941920528654516802241093365e-15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -7.65712722253446172294321996643259661334528734360491402384245881720504947974618941889e-18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.00315206806375905680927042642419402837617729582878141293654269070861386455343121457e-20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -4.67823488694203852929295147942658910943366513902587599688834805295280061153989766567e-24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.42784331030702330814913209805860341549196535898588904124741295988710093967131088558e-28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -8.48033607088050589024790469649076336393214555653483679072498587123824764496639393119e-33)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.38678975538457226036965671734240034395107611196492658679331149986658879027568426753e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -3.16374538409167657752796610509153125477921149181502343813699059314745306802450755203e-44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.81171114010755337144545373726087472870924044681995068325590922531565077520176259988e-44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -2.12518224599976943855815669111195755470966740275125423482269033756614972278324281875e-44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.34263635935093852590059444867684742593087836231966301216975070578889017560759389549e-44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -6.99309631968130547786945826997659785416033441104613988483204974783621090009511212354e-45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.97679666641473578922997037967978768575315112030381704618978104186395317308038799663e-45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.02914662367729885237175382124322264439002075945431846631322577735558916377073070681e-45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.86340999808161554465998036301386594999661959859518584115844390884231065866412899967e-46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -6.29984672974930250076096802152581521034785267762668779142172198560320648227598177519e-47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.06114476593975974108382448813738674408492134446936551295509200129676957425762855493e-47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.29100594893020385764408454293732511414412187214201378654051353911377525656418964207e-48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.01295708418471529474900073238291374279037984260859165547005675299775146592711994348e-49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -3.85645797686287223531258890899878578933576632756354490084215435591125736026369211048e-51))
      };
      T result = 0;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (k * dz + k * k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      static const T d[48] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.96888487550957841138585220622994611985377489508130630837700057227501332853802783571e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -2.54685487861380185221232939396911532213815096219473151461939579496835151047360732452e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.55552167748949946948357664110976074745481263569802725437644778906661083814390424521e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -5.96594384244972046908436642965430637119152525028640569187886872989713864303356592740e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.61206392071125071279258677015155474032235397552388756144191669459747349464120963595e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -3.26386207619544812116284495283640666091191851380376198408981809802540783149068179434e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 5.14119209063000811010153319909166529923295169312754965782286511885151365058895055453e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -6.45796781987366173516756643389460485986669456101668857494954730722833100110613651939e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.57848002216691432249235663216198359909768823927158528746723212875688277777327166892e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -5.49763614406879724077940366200859712110749040051198851331880172109188055935892153627e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.79878440118422879619706786668472072947762184485364412917440713970863702012601831523e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -2.18117881690802831493078756446312359535362496025568649851035022806893921424585651770e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.04346085879989815616410913816129789693536068884338152664460700712557690664496356491e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -4.16203393160809415901576300558568214071599197484753461939177893151363492512491744816e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.38265924315387281787206504282410481919445750548305464376847792250029309946303762942e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -3.81511194797468462062874900113250266465621301462447801385359808929965375229616323251e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 8.70467183357417683453667167003713608362889484387619452881332522389524770233094455619e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.63213714284756057716179081387804854162483375452448157161319547511874432644723929936e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 2.49457336262308128548970424731422625668341066112088155314174120975542287189333019758e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -3.07620269531262256257667392994057559249714693275375776469485257399496930680150472549e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.02187366435225062095438319574039637826727962739349882824178780109324022923104253736e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -2.32787467733840531952793152972131223922164668925536951733638661285196165517774652382e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.37936454861394636073464067115999570083944553260438583828963060925061101644510914088e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -6.13900038094091115267314860504732881894801020942472019504308090210210731167345398552e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.99245982152457111301991970385956762536902876746714853773765158096202639164545608854e-06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -4.54442771789078363056076138819579948120113702520046296917694187019851573438477912893e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 6.94926914964866053119681538255564252500374927699409988373657464226793796472537754598e-10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -6.70422350226445404683228568613754929042955139666562491583132456753472072839819516057e-12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 3.76497022572067913016091674317545685518824311681605875941103548614928596445142654530e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.10315496743460248468022132288915711635897058766184880707096827157196752704473690191e-16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.44523155331686111418148572564943890979552960547133940740999458361025242160963413278e-19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -6.73988808644577660977639361852731370450420349068670037458430193432726909116710442573e-23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.81984174508721338476616177429473006273759399181666460057479840823572355207273441830e-27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.22175387586288279410299617080503779059482571464350109824108120007426231737351306256e-31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 7.76069628266069128881232691561697202264350068700303406802593331217743381427492776348e-38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -4.55797760012119196385166534789617212222293557097685410066659856245708610260816713241e-43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.05080524465178727911313037718933906722430509738224888801936471543177490693203242080e-43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -3.06172965819221083679631791380848228306721487317680820361950045678791125133333041491e-43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.93432331242636842233461846275246262853067820473357196595254072710231999304513059821e-43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.00748867278864769255014599872258605211683566333279390476878762044733272684660155341e-43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.28864237743622222407442194745385245322813091389014701339281732300960497092852607932e-44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.48268166001867959355256619834250603827306390609511540411093653502476081933431007754e-44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 4.12528729298048868682661828587192615782099352227112063181739468363374124346195882740e-45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -9.07612869947754202550701651568771230782100517153410159160000452433811651633069544830e-46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.52878107636588388950021856856509972417665723455977954669464087724444519346925778688e-46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -1.85993987583058920731606562273777225552135287406894580374509213860879421694626894310e-47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, 1.45935754590553989061937536984206183886886144778541444482827659230313539351404405283e-48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 267, -5.55596198187619704372842019685051223079581761704170979479661609332351755409245768399e-50)),
      };
      T result = 0;
      T z = dz + 2;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (z + k * z + k * k - 1);
      }
      return result;
   }

   static double g() { return 3.71521093750000019895196601282805204391479492187500000000000000000000000000000000000e+01; }
};
//
// Lanczos Coefficients for N=58 G=4.7377167968750001136868377216160297393798828125000000000000000000000000000000000000000000000000000000000e+01
// Max experimental error (with 100-digit precision arithmetic) 89eps
// Generated with compiler: Microsoft Visual C++ version 14.2 on Win32 at Oct 14 2019
// Type precision was 334 bits or 103 max_digits10
//
struct lanczos58MP : public mpl::int_<334>
{
   template <class T>
   static T lanczos_sum(const T& z)
   {
      static const T num[58] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1113345018281139401900617648632565892997531886271137869532357720305909403460902551816191757935253298297e+96)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3306998523015859302296206498309765162943025964850544682652937543993544826195909937055933931442768067675e+96)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.8268352474100341852641008414682924186512501342670048541149232501353823124314085480609221821688640650080e+95)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.0141386651493209831831274089811845202892902987250465018037980182315440505121630111661568734223655786214e+95)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.5471258592914542867143021588167696855338160407933457150021303941501653183920330822783098409082225352649e+94)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.9029900619671536896112969373037726527696919271039767737722503312395772965346303367519600397313792277727e+94)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.4640693221647807642599419508017238730661097576373283904791677661939143319337033201070793440010339999715e+93)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.3008318491825889746629471128738097186407559167000216981000503166575445677111681454415332001004453692360e+92)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.9581961997681479620796339921050983411541681934353890086635121221338490738854914613804353855734829750740e+91)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.9562812093757389383361734216065083810639336069848105427428961256383580773444097085298844542948261561976e+90)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.0204448099875032766488461246324914875195022845573308192494619777200657405252573402181380342443147295173e+89)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.1967742477741068422882850261828584111953166907106815923693978537447964053382494347755201222673169434806e+88)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.7934308426826323386245296144249798411531042807557571947289014719343462128587532573403065724298496762577e+87)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.2112731316363645590406244691358512085700135781552933145962419229534561671703452680159817420544888863655e+86)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.7792875367286636291574914918267428357586601357931867786906344130539973047706286380462771311772443260655e+85)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6729889167781214645382609228554520670924462597343847547527634847733870633193570419233293354907537898968e+84)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.2212896890317209445182351532165347688163702283356828888187120214401748895306876396481747000574259149233e+82)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.6696330100177872775181285832060367195698947756070274789550514016487770556676430692662838020948055174819e+81)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.1787790436987713717592679695802576546659292017904142826938775295267502136179724011234775858106585474335e+80)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.3897471960941486230091766140467565375262975446763006777927562006563131703283610364746875620601313483057e+78)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.7456216276193058883301725852613487104982510349671722342387971953607146505494372335184541893296983766398e+77)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3855091017101453942895165907767739025127711826384743418432444735139729009897407793448891996393706658957e+76)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.7596919543541505311176214833582242571197674817260097288825507658966682365463135538963387631678553863542e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.5205308568827525377048378465545277453459194177597465928758908233213456983481011601239823019085100837607e+73)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.5219482007379940303673464060102227421305538751659464891709211634879751289944728765835405235116911964067e+71)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2529929196659743450833612479027982672526994601606551595908472901860777472333846227605233498599036089780e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.2371343799933848311703711400651007221131224388690941908721236718964833912418768528725638940635396117593e+68)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.8015422120873735048084230762717970714811562909437854684093754471976926045882197320273930142628411250042e+66)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.7544934847649309554128203553793825206530899558440305080175611169946112759760674627144360103537889679591e+65)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.6825279122107455550724823529764700516671135432239079326787536307050166777000218851702259835300111357922e+63)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.2137862491643897667607321980849893924427698956617156800113142641018717095876369593791167006646811605785e+61)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3186571319674830310069793304531422675865920550412368667057521060523822178148626902742473745717612071248e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.2485774541577495056605563690047646700517885331946681444187667959802672435192052254156475697343336771638e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.5749708448011926909660539105173117507305781970952738170400845272354829441805660173686553510885012106315e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.2957728089076810011933348982855845071484336995217889625353037764527196874992238794680831083872659220534e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.3029702073166140154775964636478493296495204191370975362621002042260278819606865308080233433968063616502e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.3651754658248543289835755041517884667332033286508690963566578707866450840862377843752246373050313956333e+50)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1153610110919853077136756888034644015689784357746290647672760493336488593590043112457135368750990637667e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2317697361001624160827173449861094521356739436250957059660912760537330301205792128544135654641326400704e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2591336713897064687825047779362333153357928528455245343833435858589598201567338554890557167949204480368e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1888395650074118534283640535831622534704266864697448915278695723452804398372734472806777342045155424277e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0342195380047963732482565081644533382537722066923561121178973781344452277931733870503748282543943428982e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.2659603644818244197377849119572544438113712679642935412688508428331533577481416745415954979034017913781e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.0493876469565505074955014266411916314577665137838261525037575218110760155069631757821913430706010515995e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.0380102035661800729809170954585443393379101722448652003355183633039922242269533846456702077836461601594e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4471701800736812643301636372332435764064729819100848382778186520311589961354650397132342954205946278402e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3391793076015860668513073429148720023215553662036803736563341850799950224212365984411047979319730267348e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.5746079755559823667940853958792940224952894096446926758424644873786146686411034227207764327098030858365e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.8730952802737845263408235958198085714879899902329571556827316641606659145317043632556281422581317630775e+25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1068860106460271425431255001071933519142520498723285368106201466148952145283909421951418316204387374823e+23)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.7146191714834483705600971071352976348572227109511068140706241351872914486137201205172098597005503434574e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0693465808661148698636384310734985694414606785941298207428612238748657882890467389962396214120403612061e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.5877812671795420049920049903121324610843253367633375770943367999226546473657225178995404939996916553769e+15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.1199728232559097639450648660234135896472362777482676055586950665599850474720777665093211412423570836354e+12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.9536149625531631518225658997937689548178270704423376010017262234330926444650292977902802706698285898435e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.0978434212942028068922584596812569730178537381770594492536210144311893565241511848298648224631091914130e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.8136582389003336880560289199735195790527296500214034301547642961678068653936021531092803113827675512051e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.5066282746310005024157652848110452530069867406099383166299235763422936546078420094024669537723756981791e+00))
      };
      static const T denom[58] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.1099858780486345185404564746372494973649797888116845868744704000000000000000000000000000000000000000000e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.2787481989562818154790285585287083447455750564964755439725051904000000000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.9814423801545723353305894697802330693700024131737540066597268357120000000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.2379872707088650213961476468088967648845485386801302969710798700544000000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.6028610849707118786878466090726746922395786229117545420266938564608000000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.0438085338641386367983933867412545578566598680906745465073101805977600000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.3492380371843850750819852073737242436575952022390816483000346566721536000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.5102990983560640009327873678307954320470630605324516076646299828695859200000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.6706536806860658574451791386582141941677796433514743219241257405513728000000000000000000000000000000000e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.8043159910377823217375261881225863717725312218178043012454139087028224000000000000000000000000000000000e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.9335160303959039103395501372789994167981664930485977399292586960638443520000000000000000000000000000000e+73)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1723691196565373582067598608574778201485664959896330847068435199126667264000000000000000000000000000000e+73)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4438352206893499975371296011596621640882812103203541286914698812249866240000000000000000000000000000000e+72)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.5035118941062978428961903224471555396129604077556196563177106772688896000000000000000000000000000000000e+71)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.3850308579063636618912025850466431356280013115826126911230283854233600000000000000000000000000000000000e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0836900819838640440767313247219165663768918429617679815795491004563456000000000000000000000000000000000e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.4298706038662999363673291534756645260294395901949985216872998030991360000000000000000000000000000000000e+69)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.7034492632938290783198763770876160868300549275946097010473332643840000000000000000000000000000000000000e+68)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.8389252775111283069889631347879462661235829004870041197088140774400000000000000000000000000000000000000e+67)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.8045026245181087909282578887424942187415211828585867350357833779200000000000000000000000000000000000000e+66)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6139415015963261948795856116734033680672939822172300204078108032000000000000000000000000000000000000000e+65)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3187949670862787149323475106066210874520928766305564723661651200000000000000000000000000000000000000000e+64)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.8652947935638867544975258598123447770317269214048194873129760000000000000000000000000000000000000000000e+62)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.7677767608518504817918794029244339508807996285424723220104640000000000000000000000000000000000000000000e+61)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.2641467556495415387682543812562093348191054073538914270730400000000000000000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4706517237142185424284275343767027578875512836847707758280000000000000000000000000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3177422716201009039921461100044177449498643771855331909700000000000000000000000000000000000000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.4750611894960504179018603471357527802382408288670161342000000000000000000000000000000000000000000000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.9330478552107545132079399426533304480054308936780642545000000000000000000000000000000000000000000000000e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2252899213008209896248783169149861832238011945022246000000000000000000000000000000000000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.7217402128877734970019794026056267173695376012268710000000000000000000000000000000000000000000000000000e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6785124525594845202965146591601984920001061237741320000000000000000000000000000000000000000000000000000e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.5034829951836538559219701088754272365603191632587000000000000000000000000000000000000000000000000000000e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6637415534451225329152181513653760874075278100800000000000000000000000000000000000000000000000000000000e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.6347190324222080575971451540788329488139570600000000000000000000000000000000000000000000000000000000000e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1887997793364622891010571019581486512880887480000000000000000000000000000000000000000000000000000000000e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.8047716936655096027980248782543659943937913000000000000000000000000000000000000000000000000000000000000e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.0790225058170194637643585929666817947860000000000000000000000000000000000000000000000000000000000000000e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2084757755139694726502828251712490452850000000000000000000000000000000000000000000000000000000000000000e+40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.1993406165887949960123200495814460898000000000000000000000000000000000000000000000000000000000000000000e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.6561413270081694853593737727970949550000000000000000000000000000000000000000000000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.5369993949587931998758036774526800000000000000000000000000000000000000000000000000000000000000000000000e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.6151991880839964140266135934090000000000000000000000000000000000000000000000000000000000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.4760566491969897737280290058600000000000000000000000000000000000000000000000000000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0621766041618690022809397383500000000000000000000000000000000000000000000000000000000000000000000000000e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0668475851162648374876202000000000000000000000000000000000000000000000000000000000000000000000000000000e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.5410679744103301927030500000000000000000000000000000000000000000000000000000000000000000000000000000000e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.5395191798636166637180000000000000000000000000000000000000000000000000000000000000000000000000000000000e+22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.2147541393617700530500000000000000000000000000000000000000000000000000000000000000000000000000000000000e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.1197409850631024000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.5899910954161520000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.7660425636804000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.3384219373900000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.3043596000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2435500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.5960000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }

   template <class T>
   static T lanczos_sum_expG_scaled(const T& z)
   {
      static const T num[58] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.9525834451195623352844833708766019251484824420589901934439764663785971607902920873854007890984226025387e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.5353913226536305547742063023255223178577430351061832712694094731916543407128432857596013878418835545541e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.0794265040063864629969189140227704423008747481150565046014672784314863153912417986968014303182196879938e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.0079363228400099355108372825384509227179448653556047031855632715439781114354847693263728085548142549876e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.2707926617939589278984760255796767920967443801395440059906534121478842405887081438764736244178311493742e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.0558467715603214792617326562510713807234351435955730720590956416562782623810902927804516780978093617536e+73)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.2033080198136382677487341183423922900647133007317553822084317473855699622210307747112227636535549617309e+72)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.4083202075984620999390260696805188384092629914880295342492703293697714617999890327029040466923272486691e+72)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.8486472680848031160157166025706272700120860759432294711200940780548927766558565495831077710008357816723e+71)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.1138175900123376987163744326763008390621758864911227853169463977789770550938948181501331329408292888705e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.1308645173446504769606114641352353002279158988307024846586220931944095152963117782251173128553590002677e+69)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.9120324679281068850671805700697238886511656051607969710283305054986370771314011765262161774820057407221e+68)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.5391934623114467787956807868296099285488891013887607698670475417480968323895750297127042280725645864237e+67)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1188472337439846560727348944974065682611782369759313015986429928637543111162293072546541343053710871237e+66)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.3839859706266764054133980098731552638450995195061568584278432933554090666516712057667084959928750476621e+64)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.4447818108967399523602631725512074937872788923328188245441957349314890139383300142665492699548078016636e+63)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4499038978543131336483547267337201487870950383856380855972848092343622627255235050078874464565496828787e+62)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2406238713440761680092008034709525627435037404776653076097934320335762078719054009676170121443275896208e+61)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.7885604418978046593206212118683810887087891637289025925672325034976573417819599577035132403057170493687e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4946595358499347633866345682061171079599164941304972231143049023731630101503606670560705515260702923050e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.9513336364509343072555606544850997149882112871146960235297930573184577313990835201997736426755771031690e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.6810080403717774179828377811768440096006474456859136430725963485797286723671499145937476199781240726843e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2645506501577530162354586527854723438171436013868785795302304165655061134281691393986730787208175906833e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.0397326173536285938710693563469993609385607650335155744195392582126325319144559249104863671555936627244e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2013871048927633852181093624593869925193443374649283179426834256685646298857112848793626713837584764108e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.3289402473981501234144402953028274630427307351359919167469712295902671507115725494100349633541940933204e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.6003893195733192859177412982185685700998370832033861854329854490539688708839874262744736474012759445436e+47)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.0727066732760648460132569879963196679351152547539542083857690449447906448979182035004975167951264229374e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.6613224093786683644702502032231164856098499632779487718555251290733786904948988026020403761370171127380e+44)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.7837068244515233545220503730640329863081890096623187288769939552653286719974960799952261791747246686408e+42)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.9165522010589216370481372290739111949387650621903781998282461665273890471374114621564713215060999778073e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.5033963322756782307786713227898388231893223357237704976668993902202940789294809098147971301691772648755e+39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.9740002270190565366640741207315377496247453294707036602010614289735470303392852115947439612907179897898e+37)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.4979501813196320177935632677278886186231320810970746337833403723731191512041559167496248429116298871699e+35)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.4069761263574572480694806529063097871169957005541147202721661005052475351188525148368053997631136759790e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.9402465143351975112300028774197323450509816779255592414045815705522691518043705429706146342655130265366e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4881313407932962145111410344090747527354995727266179099775181263231061106933689865989301929428200117871e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.9632810384855297696182113289305607179489664745076688933167066674562917487536504070129264115660678482556e+28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.2725546854039254943463934183457508258081242388204842730157719983831445104894295952204798603845536120509e+26)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.3452549409942328202733488146109233160716224388189413469297783823577482916724384652184656289642575364301e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.1584981954307469031764309294640361429335538748027906889408285804565757617154922327898215115914925102238e+22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.7477051072463293412749200465569024873771231987242784210248733085001441979676968452992807610739566040156e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.1960928676322401992340089643560231337218499440045721085448636142070840811937842162960498037401931828261e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6071958343894951054088987834710923196788323800906338042869660264064639861175574231931711270452905023878e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0728148958446892912836806181637101367297983133466225777574815613761623459000110314136961372678097371209e+14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.5016195836538037540865015222440859129132562974585463183411281203357391116994585543064375504964287033151e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.5579194627421707093560391662678700402228252243270030776955569538161088100804816540797427121002521848777e+09)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.7467358958842104945974935265631802960291549549388120971605962335321138414180962067591937089645593587253e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.6332135345687020156244092536439323519014928693505240553789259043703120308196994214855462894182573058664e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.9407647340128848401194943791119430377408353392667765444895085026094746239619902099149777694224783951328e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.8689665916105177080976558363071196803124010350312199134246436579120098844674440713135687985135098450472e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.8410303167649082918635413442316246381709300383012975083202578569099476956148029095902566747899894260699e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.8751938471236167190844006303614257144787097468775849169694452270347666880072332030058767126104666814717e-06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3602697452963248452045411104106027969968697679097453973063118187262717386356657326020477164267694748313e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.1131092239698914177951080885440635365109643513308304020521218234027520232780450142616357759194523396831e-11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4171068051802330111071549877540818528598812224565937952488425520621143590937182779724061728889603469712e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.8102465534711886967772860840717324393645550944967932350301450628070653877586909697550907928519183677886e-17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.6595873110837641575235810188999693348770063449340247652053167124146826554900758076335155425689417108046e-21))
      };
      static const T denom[58] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.1099858780486345185404564746372494973649797888116845868744704000000000000000000000000000000000000000000e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.2787481989562818154790285585287083447455750564964755439725051904000000000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.9814423801545723353305894697802330693700024131737540066597268357120000000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.2379872707088650213961476468088967648845485386801302969710798700544000000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.6028610849707118786878466090726746922395786229117545420266938564608000000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.0438085338641386367983933867412545578566598680906745465073101805977600000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.3492380371843850750819852073737242436575952022390816483000346566721536000000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.5102990983560640009327873678307954320470630605324516076646299828695859200000000000000000000000000000000e+75)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.6706536806860658574451791386582141941677796433514743219241257405513728000000000000000000000000000000000e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.8043159910377823217375261881225863717725312218178043012454139087028224000000000000000000000000000000000e+74)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.9335160303959039103395501372789994167981664930485977399292586960638443520000000000000000000000000000000e+73)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1723691196565373582067598608574778201485664959896330847068435199126667264000000000000000000000000000000e+73)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4438352206893499975371296011596621640882812103203541286914698812249866240000000000000000000000000000000e+72)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.5035118941062978428961903224471555396129604077556196563177106772688896000000000000000000000000000000000e+71)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.3850308579063636618912025850466431356280013115826126911230283854233600000000000000000000000000000000000e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0836900819838640440767313247219165663768918429617679815795491004563456000000000000000000000000000000000e+70)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.4298706038662999363673291534756645260294395901949985216872998030991360000000000000000000000000000000000e+69)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.7034492632938290783198763770876160868300549275946097010473332643840000000000000000000000000000000000000e+68)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.8389252775111283069889631347879462661235829004870041197088140774400000000000000000000000000000000000000e+67)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.8045026245181087909282578887424942187415211828585867350357833779200000000000000000000000000000000000000e+66)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6139415015963261948795856116734033680672939822172300204078108032000000000000000000000000000000000000000e+65)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3187949670862787149323475106066210874520928766305564723661651200000000000000000000000000000000000000000e+64)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.8652947935638867544975258598123447770317269214048194873129760000000000000000000000000000000000000000000e+62)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.7677767608518504817918794029244339508807996285424723220104640000000000000000000000000000000000000000000e+61)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.2641467556495415387682543812562093348191054073538914270730400000000000000000000000000000000000000000000e+60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4706517237142185424284275343767027578875512836847707758280000000000000000000000000000000000000000000000e+59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3177422716201009039921461100044177449498643771855331909700000000000000000000000000000000000000000000000e+58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.4750611894960504179018603471357527802382408288670161342000000000000000000000000000000000000000000000000e+56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.9330478552107545132079399426533304480054308936780642545000000000000000000000000000000000000000000000000e+55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2252899213008209896248783169149861832238011945022246000000000000000000000000000000000000000000000000000e+54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.7217402128877734970019794026056267173695376012268710000000000000000000000000000000000000000000000000000e+52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6785124525594845202965146591601984920001061237741320000000000000000000000000000000000000000000000000000e+51)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.5034829951836538559219701088754272365603191632587000000000000000000000000000000000000000000000000000000e+49)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6637415534451225329152181513653760874075278100800000000000000000000000000000000000000000000000000000000e+48)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.6347190324222080575971451540788329488139570600000000000000000000000000000000000000000000000000000000000e+46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1887997793364622891010571019581486512880887480000000000000000000000000000000000000000000000000000000000e+45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.8047716936655096027980248782543659943937913000000000000000000000000000000000000000000000000000000000000e+43)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.0790225058170194637643585929666817947860000000000000000000000000000000000000000000000000000000000000000e+41)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2084757755139694726502828251712490452850000000000000000000000000000000000000000000000000000000000000000e+40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.1993406165887949960123200495814460898000000000000000000000000000000000000000000000000000000000000000000e+38)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.6561413270081694853593737727970949550000000000000000000000000000000000000000000000000000000000000000000e+36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.5369993949587931998758036774526800000000000000000000000000000000000000000000000000000000000000000000000e+34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.6151991880839964140266135934090000000000000000000000000000000000000000000000000000000000000000000000000e+32)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.4760566491969897737280290058600000000000000000000000000000000000000000000000000000000000000000000000000e+30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0621766041618690022809397383500000000000000000000000000000000000000000000000000000000000000000000000000e+29)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0668475851162648374876202000000000000000000000000000000000000000000000000000000000000000000000000000000e+27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.5410679744103301927030500000000000000000000000000000000000000000000000000000000000000000000000000000000e+24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.5395191798636166637180000000000000000000000000000000000000000000000000000000000000000000000000000000000e+22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.2147541393617700530500000000000000000000000000000000000000000000000000000000000000000000000000000000000e+20)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.1197409850631024000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.5899910954161520000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.7660425636804000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.3384219373900000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.3043596000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2435500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.5960000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e+00))
      };
      return boost::math::tools::evaluate_rational(num, denom, z);
   }


   template<class T>
   static T lanczos_sum_near_1(const T& dz)
   {
      static const T d[57] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.7428115435862376558997260696640294500966442454383057764865555979364261022425976429085092772284462056620e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -2.9098470636910381824837046571049540905024165376037355421260734899535717859205986796003654941175663049821e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.3238837075893790405047962895917439384623673303456004530561387555917370190498457070337183015927777421697e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.1818710867341557832686340441113771687173935259473061186027784614571488561024113573068868892276452685792e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.2993023860235266500362186023069966658235478225844203105044378872946562122380318851444839578192967323002e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.1912180368478885109223112981184136172284705131966546680563227042436992258659147484337283480302880786198e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.6139396362168514560020167846095021752825297585055554665973631410697906644150591916500251626224775776603e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -4.6634802529570780684302635411632208905442609915247254971321402103134398787017804507333676554065934886041e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.8912442659385178766302183254171219505160153159736231877672680444314286802452078408167890642564966971613e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -8.5495375341789279888846473669676987489917192922529063440714880612919546655027506056499234132070763181503e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.9953131705119738082952990229447045802828349018683164227185193391425402046621590400831377644280156201943e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -8.0869264878440576629536295016426099545140305890893332059172420721937348112425650207439043097833689822279e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.2468943229845282686806294869840458047361754209121185638286960316366706592206130567941492816463381946909e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -4.1630054452110204332042164710368677551444071577885528134258060122924060097997330723219233292551334541916e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.3999790605522096466954112874401964638803242076322683502391870253501264993933977365040197484368089793844e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.1989159366950562464845132094716096217280130489202748985857413544679377185401912495211215199827572320658e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.1934569461829331860249338515130842345952840084924773082558044904616988546739452572774668343010410556738e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.9504509773595632752455387440694137287105978110568834596383152715430739170438241704706104769027458757024e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.3443645806133388465801728173699557032879027729815842645734341879259895367471331502957361274245333765268e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.7841389657268783911632604958857419096159175608100099391870386517930613345705029477317961997797817456392e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.3263071690820256271926147120746615152585668319664536987434389730938667811945831461118672855948398156120e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -9.0147093132626919930705734140711481562906727339902047165060887486522943229634999363817715017657681660287e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6071102334369619355078667017356321869894503980205853863658685021270881486705453945159213607404060788769e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -2.4383630842681764529222336922814378927628353308340684972787496418846394508904604721467821529575751627808e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.1286133930578257063116405165630725640730704580877093703762851704304775232899279568326548372516108061737e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -3.3692644530387989631839105618463769604011314569959807313957135524319268073994325844277072212634894718838e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.0185001434002208420172349105229077528666279911019323889220715333407026262532945659702898548882557288755e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -2.2263062380436433070464690537097677591705999260831413712451347815254366329816629201276858957857633159301e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3353139137616007516399230327860238131259262384241927446052286747295202825849286263834250215682608257253e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -6.4195296490132509971486202730620553572364290374630355485943647261186824263045109480840776574102283117802e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4317122347228636250935502644177076353230092364254437169086888072416523402594888598203831115818706298660e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -7.1116401895712757982110113602984617019673481361481571387831251493876840058692034827601429496435585930404e-10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.5670513725610592064219547014678059161831695745208975662601768645833855387548399162792665797683105415471e-11)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -2.5260646573782426214538742421235036428884108667170155792319000553615370065204316405736286181192562052741e-13)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.8733941140402994491303889414468143535352703737984260742434438594748723811519620261843045417069155536008e-15)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -2.2056577386533803436177435861654503857156683256162837031477537088773467287029340043192251326622052870465e-17)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0800931908111468667692114269177968638854301916581458483888602672737990097921789711610932672475861001994e-19)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -3.1394496063229749361320583264730160299755358905685423988290847461006304721866926670065937080941317261860e-22)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.9274792450690899704710418936228583371317511924974285800543619049655797337986789330520046119470213637029e-25)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -3.6755223694460581147198804301148887855673388834174147952660909117310772380696272509534022412162907695615e-28)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.0900225572290190698649005514102451202434886376168031607266730359956445823401555519988657304185446296554e-31)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -9.9024654365441546862123770751558932020185466144030552444525982172062830502919077739315180881006251414076e-36)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.8401409342622593684683418781548298452423138616374012775158698166549875044639516398318322746270349752685e-40)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -3.5509630598854356990909041057867643490325454268799807828321700754551551638916414252086883229952366535870e-46)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.9657368621557895021941755410664403576678059155205990324521617023894341207730153756506985102036313160050e-53)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.1241524357510805548447946933511907560696760445486904751298400852883111381333879290886698299626153350908e-56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -4.7964162808880362344077245029516705444495147547810713281233063336303852616170950226238109362470235010501e-56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4347533225708328774658804669618167169988215409141087123873615495103598105335608364600471221104934664685e-56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.0479456512751620605515259329413452300868549663188231386090986433707412031359388913176316853836715844262e-56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.7702413521297716432095744498633554542922110264197418721766960511801656053258128107225421844535437414965e-57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.1172472608352701420007786830697096251191397605488453881854169564652213625087310098345004527869954910248e-57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.6812745632569367999094411216914045337299158629329933113885292569194336792226670477775358040058394820261e-58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -5.0940899331938988550733986571553993560122282178667089274000678435292023664866549561623149093538369388095e-59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.3990413208296530427672464741601390159299766279095164840741566288983363882627013060272921361169788365065e-60)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -7.7457129143553309783200493485805274806723421714099516756842649457968744562118823820376006640033350414786e-61)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.2197041624023821332095587670691140585500743194040476924703015039327636848179534186153838408314119752602e-62)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.7047723462914761686649775942852334082321432815609988664912617018852671864913042227323649381537503179251e-63))
      };
      T result = 0;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (k * dz + k * k);
      }
      return result;
   }

   template<class T>
   static T lanczos_sum_near_2(const T& dz)
   {
      static const T d[57] = {
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.1662905928411697905258405605703584854097058546117590807018322943510409180314506715189727704455120038255e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -5.2865276330524469503633315971419666291446796012439883450060362841602955272110829770969508625348151344110e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.2219660233923727011658095806761811032946810239961697070082822527681712138952842835075740681056558805545e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -2.1471898769829091933140257345448028662421485176900416920959625594396387278224755880120811584841740984425e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.8108506629663017023576288202278104175560258464862168051112507988271579950339846330515110250714726935097e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -2.1641711509053717182383657533380473756042975817599106212749608660372803916254884680397716414106410031258e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.7489314096332467747675642425312337606864580547949737826112168114861342185432170188369658631861013639815e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -8.4724786849037234839356479194796060821229710729537665905710787287944594232273268923957076043820505152123e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.2519817172723530341690766298268651248835138906207131843033724432774951682486028564360090931772389689559e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.5532557359535150161664447223551521375257233550726531311510295929236015954526292484349276082718903034254e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.6342429895113330905046033016646335037353665731613082814708242545769972330151611687916571207590466290675e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.4692098728454232462170975109580416322729103219472916098312631671127459727620998375190829900768135597920e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1349180467693002337442162064332669746019375174297711955789632954762786417078957693464106665366428215075e+07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -7.5632302457640823075347957970491681679223652407676639065204661139779189254891918422939089928834865040327e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.3602139028787279699176338906645576475531513654199858202869033753097108161316098377414795783646795195944e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -2.1781564770643699404196555460989570673904909639666069579972319801392840126860473187109481731392038690946e+06)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.4353253130211266653540571778153094858600293629830279655688388322741879148986220244682882658942335696445e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -3.5435240282512309592551076592408102377879926851607713336111442642626445924365822019934119240733876795141e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.1526261667865026060062626782458367235346781022115660163078870751778821943822555542992468135865376114967e+05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -3.2413730814968410697200758155304075536237935327710818457186113369696959175294570307981528264101370866348e+04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 7.8599121870734324248863589686279377309239133065805881581800599402339736200679855503773749459170204132791e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.6377668257261558531404127005983837331367557885561101077531977908963614729639690754651655156526373276599e+03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.9197523005380738606902261922935494412976764483416941348495962454277708460181544072006267228271644673220e+02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -4.4299489087401028354837196175165054006400597113077360830473379437347868849689360216846724154620950362855e+01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.6839760968600182694175565648725095464923929298411152920412000388870216851203160830044900161240884736073e+00)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -6.1211841186792218837353890586540028399627732652134291573554656108241153104326807995702108419307165078338e-01)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.4839254672775355063516929562453773790762108550565402409275434242243168512055903461692839292134763272349e-02)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -4.0446900436498039824930882038167913752019117688189621941029942455212255095042589886820040572146790228293e-03)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.4259604540679199219674266589960434988970631237379997092121152488502215342212782655207033499356597207732e-04)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.1662819432736824144566608258931608012907080530838741239388324423687927670373421600053389501450361346193e-05)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.4178658338791266451928074758746505934456557516996342163398372684287698861436808924171077524254379083060e-07)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.2920226237184384323812413019935857990607979507557513737013812388314162555859654820969350110668735340875e-08)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 2.8469744980165771885329432448423013838689752067078815649758371478630140330027362965330275099502650356588e-10)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -4.5892826398814320794517525606378651211265264604594412375068674574147845337527175224917342725746156587979e-12)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 5.2203009477948187935105053137418678388241808213010074775815651125432941229706488789085110080538814893051e-14)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -4.0071764354709508994124453699849221531650525550222671429443364906737525324879927767883427021892836113984e-16)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.9622826817063220108064709012206725989960307918188611950534485874633362779114658729325037773191270734438e-18)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -5.7036630218461017970632409261270377964650553992557669797810281540623485979537248227679839983488911988094e-21)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 8.9521045677595146905259180949264876791866636312953998390445326640202924376790113203607493358448772833037e-24)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -6.6775848168912941217401852995810450906306734023598389296925437100201008270660133423441040356945585338087e-27)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.9803220730550206978978787593484142931342246644534949059896766020244771367805337770466949634458076040872e-30)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.7990518408632014329157379404081596938445109431006332980973114098434711412087164821902160534448534060827e-34)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.3431158699278209062613044372715762567352361096055698830148755308547762997581572887297656338501121085689e-39)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -6.4512889953126515222750706142725047033542840408835574821053782432341231211557553861105639807627133237603e-45)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 3.5712949902991134825521996829159928713052863755852691744249239606303690103524531369836280748984998481830e-52)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.4759729775024569614177814940895080213705811461280570147704308137991423094447188671055911764625031287276e-54)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -8.7139930908855290554513295892599839263800779861752158476251019536427498588557825717368115652274945883813e-55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.4233907960475632932218656634142149345608818585031048657672084471779614927408553135306513446078671755098e-55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.9038779434600254520302015223066920881918842711110146920259738288262023701119647715089263879918197432213e-55)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 6.8496675787588211884718935004472008052724911912938728434253127344504858415646126692300012586778132949530e-56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -2.0297831425771921486916100519977404266947748691416922167961096320991761190057145636867796238254884443749e-56)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 4.8712635957157169685941764659378666193084431023269775678410920776942423498574161990151090862987707496709e-57)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -9.2547981414954593575052648297803007916203550348517241024188916829704035005537638171060497453987203484125e-58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 1.3442368462845102652888331895188311507347132911248615838604886657492687976498041391593065929847726694329e-58)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -1.4072191583666827552083642053400837423373421659342004966749178938352282349787956685790064733145823615594e-59)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, 9.4830105111768027064991929773116440103527225991506525293838437779387018237471847818347678097555938325901e-61)),
         static_cast<T>(BOOST_MATH_BIG_CONSTANT(T, 334, -3.0971820578438671326698299395846141298596444280876893045395714232857615676238866889122426000206260012404e-62)),
      };
      T result = 0;
      T z = dz + 2;
      for (unsigned k = 1; k <= sizeof(d) / sizeof(d[0]); ++k)
      {
         result += (-d[k - 1] * dz) / (z + k * z + k * k - 1);
      }
      return result;
   }

   static double g() { return 4.7377167968750001136868377216160297393798828125000000000000000000000000000000000000000000000000000000000e+01; }
};


//
// placeholder for no lanczos info available:
//
struct undefined_lanczos : public mpl::int_<INT_MAX - 1> { };

#if 0
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
#define BOOST_MATH_FLT_DIGITS ::std::numeric_limits<float>::digits
#define BOOST_MATH_DBL_DIGITS ::std::numeric_limits<double>::digits
#define BOOST_MATH_LDBL_DIGITS ::std::numeric_limits<long double>::digits
#else
#define BOOST_MATH_FLT_DIGITS FLT_MANT_DIG
#define BOOST_MATH_DBL_DIGITS DBL_MANT_DIG
#define BOOST_MATH_LDBL_DIGITS LDBL_MANT_DIG
#endif
#endif

typedef mpl::list<
   lanczos6m24, 
/*   lanczos6, */
   lanczos13m53, 
/*   lanczos13, */
   lanczos17m64, 
   lanczos24m113, 
   lanczos22,
   lanczos32MP, 
   lanczos35MP,
   lanczos48MP,
   lanczos49MP,
   lanczos49MP_2,
   lanczos58MP,
   undefined_lanczos> lanczos_list;

template <class Real, class Policy>
struct lanczos
{
   typedef typename mpl::if_<
      typename mpl::less_equal<
         typename policies::precision<Real, Policy>::type,
         mpl::int_<0>
      >::type,
      mpl::int_<INT_MAX - 2>,
      typename policies::precision<Real, Policy>::type
   >::type target_precision;

   typedef typename mpl::deref<typename mpl::find_if<
      lanczos_list, 
      mpl::less_equal<target_precision, mpl::_1> >::type>::type type;
};

} // namespace lanczos
} // namespace math
} // namespace boost

#if !defined(_CRAYC) && !defined(__CUDACC__) && (!defined(__GNUC__) || (__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ > 3)))
#if ((defined(_M_IX86_FP) && (_M_IX86_FP >= 2)) || defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64)) && !defined(_MANAGED)
#include <boost/math/special_functions/detail/lanczos_sse2.hpp>
#endif
#endif

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_LANCZOS




