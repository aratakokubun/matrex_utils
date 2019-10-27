/*
  CreateVec : データ列 x の鏡映変換(ハウスホルダー変換)を行い、変換後のベクトルのノルムの平方を返す

  ベクトル x とノルムが等しく、第一成分以外がゼロになるようなベクトル x' を求め、
  v = x - x' に変換する(ハウスホルダー変換)。ノルム ||x|| ( = ||x'|| ) を返り値として返す。
*/
template< class T >
T CreateVec( std::vector< T >* x )
{
  T v = std::sqrt( std::inner_product( x->begin(), x->end(), x->begin(), T() ) ); // vec のノルム

  if ( x->front() > 0 )
    x->front() += v;
  else
    x->front() -= v;

  return( std::inner_product( x->begin(), x->end(), x->begin(), T() ) );
}

/*
  Householder_ST : ハウスホルダー行列による行列 mat の相似変換

  v は鏡映変換したベクトルで、ハウスホルダー行列の要素となる。
  対象範囲は col 以降の行とする(すでに変換した行は不要)
  norm は v のノルムである。
*/
template< class T >
void Householder_ST
( SquareMatrix< T >* mat, typename SquareMatrix< T >::size_type col,
  const std::vector< T >& v, T norm )
{
  typedef typename SquareMatrix< T >::size_type size_type;

  SquareMatrix< T > buff( mat->size() );

  /* 左側からの乗算(結果をbuffへ) */

  // col行より上側の要素は変化しない
  for ( size_type r = 0 ; r < col ; ++r )
    for ( size_type c = 0 ; c < mat->size() ; ++c )
      buff[r][c] = (*mat)[r][c];

  for ( size_type r = col ; r < mat->size() ; ++r ) {
    for ( size_type c = 0 ; c < mat->size() ; ++c ) {
      for ( size_type i = 0 ; i < v.size() ; ++i )
        buff[r][c] -= v[i] * (*mat)[i + col][c];
      buff[r][c] *= 2 * v[r - col] / norm;
      buff[r][c] += (*mat)[r][c];
    }
  }

  /* 右側からの乗算(結果をmatへ) */

  for ( size_type r = 0 ; r < mat->size() ; ++r ) {
    for ( size_type c = 0 ; c < col ; ++c )
      (*mat)[r][c] = buff[r][c];
    for ( size_type c = col ; c < mat->size() ; ++c ) {
      (*mat)[r][c] = 0;
      for ( size_type i = 0 ; i < v.size() ; ++i )
        (*mat)[r][c] -= v[i] * buff[r][i + col];
      (*mat)[r][c] *= 2 * v[c - col] / norm;
      (*mat)[r][c] += buff[r][c];
    }
  }
}

/*
  Householder_TD : ハウスホルダー法による行列 r の三重対角行列への変換

  q は直交行列を得るための行列 ( H(v)の積 )
*/
template< class T >
void Householder_TD
( SquareMatrix< T >* r, SquareMatrix< T >* q )
{
  assert( r != 0 );

  typedef typename SquareMatrix< T >::size_type size_type;
  typedef typename SquareMatrix< T >::const_iterator const_iterator;

  if ( q != 0 ) {
    q->assign( r->size() );
    for ( size_type i = 0 ; i < r->size() ; ++i )
      (*q)[i][i] = 1;
  }

  if ( r->size() < 3 ) return;

  for ( size_type c = 0 ; c < r->size() - 2 ; ++c ) {
    // 対象の列ベクトルを鏡映変換
    const_iterator colIt = r->column( c, c + 1 );
    std::vector< T > v( colIt, colIt.end() );
    T norm = CreateVec( &v );

    // ハウスホルダー行列による相似変換
    Householder_ST( r, c + 1, v, norm );

    // 直交行列の計算
    if ( q == 0 ) continue;
    if ( c == 0 ) {
      for ( size_type j = 0 ; j < r->size() - 1 ; ++j ) {
        (*q)[j + 1][j + 1] = 1 - 2 * v[j] * v[j] / norm;
        for ( size_type i = j + 1 ; i < r->size() - 1 ; ++i )
          (*q)[j + 1][i + 1] = (*q)[i + 1][j + 1] = - 2 * v[j] * v[i] / norm;
      }
    } else {
      MultHouseholder( q, c + 1, v, norm );
    }
  }
}

/*
  Householder_GR : ギブンス回転による三重対角行列 r の上三角行列化

  q は直交行列( G(i,j,θ)の積 )、sz は処理範囲(行列数)をそれぞれ表す。
*/
template< class T >
void Householder_GR
( SquareMatrix< T >* r, SquareMatrix< T >* q, typename SquareMatrix< T >::size_type sz )
{
  assert( r != 0 );

  typedef typename SquareMatrix< T >::size_type size_type;

  if ( q != 0 ) {
    q->assign( r->size() );
    for ( size_type i = 0 ; i < r->size() ; ++i )
      (*q)[i][i] = 1;
  }

  if ( r->size() < 2 ) return;

  for ( size_type i = 0 ; i < sz ; ++i ) {
    // 対角成分の下の成分がゼロなら何も処理しない(処理不要)
    if ( (*r)[i + 1][i] == 0 ) continue;

    T d = std::sqrt( std::pow( (*r)[i][i], 2.0 ) + std::pow( (*r)[i + 1][i], 2.0 ) );

    T sin = (*r)[i + 1][i] / d;
    T cos = (*r)[i][i] / d;

    // ギブンス回転により r[i+1][i]をゼロにする
    for ( size_type j = 0 ; j < 2 ; ++j ) {
      T d1 = (*r)[i][i + j];
      T d2 = (*r)[i + 1][i + j];
      (*r)[i][i + j] = d1 * cos + d2 * sin;
      (*r)[i + 1][i + j] = -d1 * sin + d2 * cos;
    }
    if ( i + 2 < r->size() ) {
      d = (*r)[i + 1][i + 2];
      (*r)[i][i + 2] = d * sin;
      (*r)[i + 1][i + 2] = d * cos;
    }

    // 直交行列の計算
    if ( q == 0 ) continue;
    if ( i == 0 ) {
      (*q)[0][0] = (*q)[1][1] = cos;
      (*q)[0][1] = sin;
      (*q)[1][0] = -sin;
    } else {
      for ( size_type c = 0 ; c < r->size() ; ++c ) {
        T d1 = (*q)[i][c] * cos + (*q)[i + 1][c] * sin;
        T d2 = -(*q)[i][c] * sin + (*q)[i + 1][c] * cos;
        (*q)[i][c] = d1;
        (*q)[i + 1][c] = d2;
      }
    }
  }
  if ( q != 0 )
    q->transpose(); // 直交行列の転値
}

/*
  Householder_QR : 三重対角行列を作成して QR 変換により固有値・固有ベクトルを求める

  mat : 変換対象の対称行列
  r : 変換後の上三角行列を格納する正方行列へのポインタ
  q : 変換後の直交行列を格納する正方行列へのポインタ
  e : 収束判定のためのしきい値
  maxCnt : 収束しなかった場合の最大処理回数
*/
template< class T >
void Householder_QR
( const SymmetricMatrix< T >& mat, SquareMatrix< T >* r, SquareMatrix< T >* q, T e, unsigned int maxCnt )
{
  assert( &mat != 0 && r != 0 );
  ErrLib::NNNum( e );

  typedef typename SquareMatrix< T >::size_type size_type;

  mat.clone( r );

  unsigned int cnt = 0; // 処理回数
  Householder_TD( r, q );
  if ( q != 0 ) q->transpose();

  size_type sz = r->size() - 1;
  while ( cnt < maxCnt ) {
    SquareMatrix< T > qk;

    Householder_GR( r, &qk, sz );
    *r *= qk; // RQを計算

    if ( q != 0 ) *q *= qk;   // 直交行列の計算
    if ( std::abs( (*r)[sz][sz - 1] ) < e ) --sz;
    if ( sz == 0 ) break;
    ++cnt;
  }

  // 処理が最大処理回数を超えたら例外を投げる
  if ( cnt >= maxCnt )
    throw ExceptionExcessLimit( maxCnt );
}


/*
  Householder_EV22 : 2x2行列の固有値を計算する

  double a,b,c,d : 行列の要素(それぞれ左上・右上・左下・右下)
*/
template< class T >
T Householder_EV22( T a, T b, T c, T d )
{
  T d1 = a + d;
  T d2 = a * d - b * c;

  T a1 = ( d1 + std::sqrt( d1 * d1 - 4 * d2 ) ) / 2;
  T a2 = ( d1 - std::sqrt( d1 * d1 - 4 * d2 ) ) / 2;

  if ( std::abs( d - a1 ) < std::abs( d - a2 ) )
    return( a1 );
  else
    return( a2 );
}

/*
  Householder_DoubleShiftQR : 原点シフト付きQR変換

  mat : 変換対象の対称行列
  r : 変換後の上三角行列を格納する正方行列へのポインタ
  q : 変換後の直交行列を格納する正方行列へのポインタ
  e : 収束判定のためのしきい値
  maxCnt : 収束しなかった場合の最大処理回数
*/
template< class T >
void Householder_DoubleShiftQR
( const SymmetricMatrix< T >& mat, SquareMatrix< T >* r, SquareMatrix< T >* q, T e, unsigned int maxCnt )
{
  assert( &mat != 0 && r != 0 );
  ErrLib::NNNum( e );

  typedef typename SquareMatrix< T >::size_type size_type;

  mat.clone( r );

  unsigned int cnt = 0; // 処理回数
  Householder_TD( r, q );
  if ( q != 0 ) q->transpose();

  size_type sz = r->size() - 1;
  while ( cnt < maxCnt ) {
    SquareMatrix< T > qk;
    T u = Householder_EV22( (*r)[sz - 1][sz - 1], (*r)[sz - 1][sz], (*r)[sz][sz - 1], (*r)[sz][sz] );
    for ( size_type i = 0 ; i <= sz ; ++i )
      (*r)[i][i] -= u;

    Householder_GR( r, &qk, sz );
    *r *= qk; // RQを計算
    for ( size_type i = 0 ; i <= sz ; ++i )
      (*r)[i][i] += u;

    if ( q != 0 ) *q *= qk;   // 直交行列の計算
    if ( std::abs( (*r)[sz][sz - 1] ) < e ) --sz;
    if ( sz == 0 ) break;
    ++cnt;
  }

  // 処理が最大処理回数を超えたら例外を投げる
  if ( cnt >= maxCnt )
    throw ExceptionExcessLimit( maxCnt );
}
