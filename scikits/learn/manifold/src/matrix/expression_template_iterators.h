/**
 * \file expression_template_iterators.h
 * Defines the iterators for expression templates
 */

#ifndef EXPRESSIONTEMPLATEITERATORS
#define EXPRESSIONTEMPLATEITERATORS

namespace Matrix
{
  /// Iterator for opposition
  template<class Expression>
  class OpposeExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&> Parent;

    typedef typename Expression::const_iterator ExpressionIterator;

    /// Iterator on an expression
    ExpressionIterator expressionIterator;

  public:
    typedef OpposeExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator is the iterator that will be opposed
     */
    OpposeExpressionIterator(const ExpressionIterator& expressionIterator)
      :expressionIterator(expressionIterator)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const OpposeExpressionIterator& other)
    {
      return expressionIterator != other.expressionIterator;
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      return -*expressionIterator;
    }

    Self lineBegin()
    {
      return Self(expressionIterator.lineBegin());
    }

    Self lineEnd()
    {
      return Self(expressionIterator.lineEnd());
    }

    Self columnBegin()
    {
      return Self(expressionIterator.columnBegin());
    }

    Self columnEnd()
    {
      return Self(expressionIterator.columnEnd());
    }
  };

  /// Iterator for addition
  template<class Expression>
  class AdditionConstantExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&> Parent;

    typedef typename Expression::const_iterator ExpressionIterator;

    /// Iterator on an expression
    ExpressionIterator expressionIterator;
    /// Value to add
    typename Parent::value_type constant;

  public:
    typedef AdditionConstantExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator is the iterator that will be added
     */
    AdditionConstantExpressionIterator(const ExpressionIterator& expressionIterator, typename Parent::value_type constant)
      :expressionIterator(expressionIterator), constant(constant)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const AdditionConstantExpressionIterator& other)
    {
      return expressionIterator != other.expressionIterator;
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      return *expressionIterator + constant;
    }

    Self lineBegin()
    {
      return Self(expressionIterator.lineBegin(), constant);
    }

    Self lineEnd()
    {
      return Self(expressionIterator.lineEnd(), constant);
    }

    Self columnBegin()
    {
      return Self(expressionIterator.columnBegin(), constant);
    }

    Self columnEnd()
    {
      return Self(expressionIterator.columnEnd(), constant);
    }
  };

  /// Iterator for addition
  template<class Expression>
  class ConstantAdditionExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&> Parent;

    typedef typename Expression::const_iterator ExpressionIterator;

    /// Iterator on an expression
    ExpressionIterator expressionIterator;
    /// Value to add
    typename Parent::value_type constant;

  public:
    typedef ConstantAdditionExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator is the iterator that will be added
     */
    ConstantAdditionExpressionIterator(const ExpressionIterator& expressionIterator, typename Parent::value_type constant)
      :expressionIterator(expressionIterator), constant(constant)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const ConstantAdditionExpressionIterator& other)
    {
      return expressionIterator != other.expressionIterator;
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      return constant + *expressionIterator;
    }

    Self lineBegin()
    {
      return Self(expressionIterator.lineBegin(), constant);
    }

    Self lineEnd()
    {
      return Self(expressionIterator.lineEnd(), constant);
    }

    Self columnBegin()
    {
      return Self(expressionIterator.columnBegin(), constant);
    }

    Self columnEnd()
    {
      return Self(expressionIterator.columnEnd(), constant);
    }
  };

  /// Iterator for addition
  template<class Expression1, class Expression2>
  class AdditionExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression1::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression1::Data_Type>::type*, typename boost::add_const<typename Expression1::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression1::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression1::Data_Type>::type*, typename boost::add_const<typename Expression1::Data_Type>::type&> Parent;

    typedef typename Expression1::const_iterator ExpressionIterator1;
    typedef typename Expression2::const_iterator ExpressionIterator2;

    /// Iterator on an expression
    ExpressionIterator1 expressionIterator1;
    /// Iterator on an expression
    ExpressionIterator2 expressionIterator2;

  public:
    typedef AdditionExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator1 is the firt term
     * @param expressionIterator2 is the second term
     */
    AdditionExpressionIterator(const ExpressionIterator1& expressionIterator1, const ExpressionIterator2& expressionIterator2)
      :expressionIterator1(expressionIterator1), expressionIterator2(expressionIterator2)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator1;
      ++expressionIterator2;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const AdditionExpressionIterator& other)
    {
      return (expressionIterator1 != other.expressionIterator1) && (expressionIterator2 != other.expressionIterator2);
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      return *expressionIterator1 + *expressionIterator2;
    }

    Self lineBegin()
    {
      return Self(expressionIterator1.lineBegin(), expressionIterator2.lineBegin());
    }

    Self lineEnd()
    {
      return Self(expressionIterator1.lineEnd(), expressionIterator2.lineEnd());
    }

    Self columnBegin()
    {
      return Self(expressionIterator1.columnBegin(), expressionIterator2.columnBegin());
    }

    Self columnEnd()
    {
      return Self(expressionIterator1.columnEnd(), expressionIterator2.columnEnd());
    }
  };

  /// Iterator for substraction
  template<class Expression>
  class SubstractionConstantExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&> Parent;

    typedef typename Expression::const_iterator ExpressionIterator;

    /// Iterator on an expression
    ExpressionIterator expressionIterator;
    /// Value to add
    typename Parent::value_type constant;

  public:
    typedef SubstractionConstantExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator is the iterator that will be substracted
     */
    SubstractionConstantExpressionIterator(const ExpressionIterator& expressionIterator, typename Parent::value_type constant)
      :expressionIterator(expressionIterator), constant(constant)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const SubstractionConstantExpressionIterator& other)
    {
      return expressionIterator != other.expressionIterator;
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      return *expressionIterator - constant;
    }

    Self lineBegin()
    {
      return Self(expressionIterator.lineBegin(), constant);
    }

    Self lineEnd()
    {
      return Self(expressionIterator.lineEnd(), constant);
    }

    Self columnBegin()
    {
      return Self(expressionIterator.columnBegin(), constant);
    }

    Self columnEnd()
    {
      return Self(expressionIterator.columnEnd(), constant);
    }
  };

  /// Iterator for substraction
  template<class Expression>
  class ConstantSubstractionExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&> Parent;

    typedef typename Expression::const_iterator ExpressionIterator;

    /// Iterator on an expression
    ExpressionIterator expressionIterator;
    /// Value to add
    typename Parent::value_type constant;

  public:
    typedef ConstantSubstractionExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator is the iterator that will be substracted
     */
    ConstantSubstractionExpressionIterator(const ExpressionIterator& expressionIterator, typename Parent::value_type constant)
      :expressionIterator(expressionIterator), constant(constant)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const ConstantSubstractionExpressionIterator& other)
    {
      return expressionIterator != other.expressionIterator;
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      return constant - *expressionIterator;
    }

    Self lineBegin()
    {
      return Self(expressionIterator.lineBegin(), constant);
    }

    Self lineEnd()
    {
      return Self(expressionIterator.lineEnd(), constant);
    }

    Self columnBegin()
    {
      return Self(expressionIterator.columnBegin(), constant);
    }

    Self columnEnd()
    {
      return Self(expressionIterator.columnEnd(), constant);
    }
  };

  /// Iterator for Substraction
  template<class Expression1, class Expression2>
  class SubstractionExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression1::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression1::Data_Type>::type*, typename boost::add_const<typename Expression1::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression1::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression1::Data_Type>::type*, typename boost::add_const<typename Expression1::Data_Type>::type&> Parent;

    typedef typename Expression1::const_iterator ExpressionIterator1;
    typedef typename Expression2::const_iterator ExpressionIterator2;

    /// Iterator on an expression
    ExpressionIterator1 expressionIterator1;
    /// Iterator on an expression
    ExpressionIterator2 expressionIterator2;

  public:
    typedef SubstractionExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator1 is the firt term
     * @param expressionIterator2 is the second term
     */
    SubstractionExpressionIterator(const ExpressionIterator1& expressionIterator1, const ExpressionIterator2& expressionIterator2)
      :expressionIterator1(expressionIterator1), expressionIterator2(expressionIterator2)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator1;
      ++expressionIterator2;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const SubstractionExpressionIterator& other)
    {
      return (expressionIterator1 != other.expressionIterator1) && (expressionIterator2 != other.expressionIterator2);
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      return *expressionIterator1 - *expressionIterator2;
    }

    Self lineBegin()
    {
      return Self(expressionIterator1.lineBegin(), expressionIterator2.lineBegin());
    }

    Self lineEnd()
    {
      return Self(expressionIterator1.lineEnd(), expressionIterator2.lineEnd());
    }

    Self columnBegin()
    {
      return Self(expressionIterator1.columnBegin(), expressionIterator2.columnBegin());
    }

    Self columnEnd()
    {
      return Self(expressionIterator1.columnEnd(), expressionIterator2.columnEnd());
    }
  };

  /// Iterator for multiplication
  template<class Expression>
  class MultiplicationConstantExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&> Parent;

    typedef typename Expression::const_iterator ExpressionIterator;

    /// Iterator on an expression
    ExpressionIterator expressionIterator;
    /// Value to multiply
    typename Parent::value_type constant;

  public:
    typedef MultiplicationConstantExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator is the iterator that will be multiplied
     */
    MultiplicationConstantExpressionIterator(const ExpressionIterator& expressionIterator, typename Parent::value_type constant)
      :expressionIterator(expressionIterator), constant(constant)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const MultiplicationConstantExpressionIterator& other)
    {
      return expressionIterator != other.expressionIterator;
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      return *expressionIterator * constant;
    }

    Self lineBegin()
    {
      return Self(expressionIterator.lineBegin(), constant);
    }

    Self lineEnd()
    {
      return Self(expressionIterator.lineEnd(), constant);
    }

    Self columnBegin()
    {
      return Self(expressionIterator.columnBegin(), constant);
    }

    Self columnEnd()
    {
      return Self(expressionIterator.columnEnd(), constant);
    }
  };

  /// Iterator for multiplication
  template<class Expression>
  class ConstantMultiplicationExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression::Data_Type>::type*, typename boost::add_const<typename Expression::Data_Type>::type&> Parent;

    typedef typename Expression::const_iterator ExpressionIterator;

    /// Iterator on an expression
    ExpressionIterator expressionIterator;
    /// Value to multiply
    typename Parent::value_type constant;

  public:
    typedef ConstantMultiplicationExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator is the iterator that will be multiplied
     */
    ConstantMultiplicationExpressionIterator(const ExpressionIterator& expressionIterator, typename Parent::value_type constant)
      :expressionIterator(expressionIterator), constant(constant)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const ConstantMultiplicationExpressionIterator& other)
    {
      return expressionIterator != other.expressionIterator;
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      return constant * *expressionIterator;
    }

    Self lineBegin()
    {
      return Self(expressionIterator.lineBegin(), constant);
    }

    Self lineEnd()
    {
      return Self(expressionIterator.lineEnd(), constant);
    }

    Self columnBegin()
    {
      return Self(expressionIterator.columnBegin(), constant);
    }

    Self columnEnd()
    {
      return Self(expressionIterator.columnEnd(), constant);
    }
  };

  /// Iterator for multiplication
  template<class Expression1, class Expression2>
  class MultiplicationExpressionIterator : public std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression1::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression1::Data_Type>::type*, typename boost::add_const<typename Expression1::Data_Type>::type&>
  {
    typedef std::iterator<std::forward_iterator_tag, typename boost::add_const<typename Expression1::Data_Type>::type, ptrdiff_t, typename boost::add_const<typename Expression1::Data_Type>::type*, typename boost::add_const<typename Expression1::Data_Type>::type&> Parent;

    typedef typename Expression1::const_iterator ExpressionIterator1;
    typedef typename Expression2::const_iterator ExpressionIterator2;

    /// Iterator on an expression
    ExpressionIterator1 expressionIterator1;
    /// Iterator on an expression
    ExpressionIterator2 expressionIterator2;

  public:
    typedef MultiplicationExpressionIterator Self;

    /**
     * Constructor
     * @param expressionIterator1 is the first term
     * @param expressionIterator2 is the second term
     */
    MultiplicationExpressionIterator(const ExpressionIterator1& expressionIterator1, const ExpressionIterator2& expressionIterator2)
      :expressionIterator1(expressionIterator1), expressionIterator2(expressionIterator2)
    {
    }

    /**
     * Preincrementation
     * @return an iterator
     */
    Self& operator++()
    {
      ++expressionIterator1;
      ++expressionIterator2;
      return *this;
    }

    /**
     * Equality comparison
     * @param other is the other iterator to check
     * @result true if the 2 iterators are differents
     */
    bool operator!=(const MultiplicationExpressionIterator& other)
    {
      return (expressionIterator1 != other.expressionIterator1) && (expressionIterator2 != other.expressionIterator2);
    }

    /**
     * Dereferencement
     * @return the associated value
     */
    typename Parent::value_type operator*()
    {
      typename Expression1::const_iterator itLine = expressionIterator1.lineBegin(), itEndLine = expressionIterator1.lineEnd();
      typename Expression2::const_iterator itColumn = expressionIterator2.columnBegin(), itEndColumn = expressionIterator2.columnEnd();
      typename Expression1::Data_Type accumulation = DataTypeTraits<typename Expression1::Data_Type>::zero(*expressionIterator1 * *expressionIterator2);
      for(; itLine != itEndLine; ++itLine, ++itColumn)
      {
        accumulation += *itLine * *itColumn;
      }
      return accumulation;
    }

    Self lineBegin()
    {
      return Self(expressionIterator1.lineBegin(), expressionIterator2.lineBegin());
    }

    Self lineEnd()
    {
      return Self(expressionIterator1.lineEnd(), expressionIterator2.lineEnd());
    }

    Self columnBegin()
    {
      return Self(expressionIterator1.columnBegin(), expressionIterator2.columnBegin());
    }

    Self columnEnd()
    {
      return Self(expressionIterator1.columnEnd(), expressionIterator2.columnEnd());
    }
  };
}

#endif
