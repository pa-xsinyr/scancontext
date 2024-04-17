#include "Scancontext.h"

// namespace SC2
// {

void coreImportTest (void)//测试函数，无参数，无返回值，占位函数，表示成功导入库或模块
{
    cout << "scancontext lib is successfully imported." << endl;
} // coreImportTest


float rad2deg(float radians)//弧度到角度转换
{
    return radians * 180.0 / M_PI;
}

float deg2rad(float degrees)
{
    return degrees * M_PI / 180.0;
}


float xy2theta( const float & _x, const float & _y )//根据笛卡尔坐标计算角度，单位为度
{
    if ( _x >= 0 & _y >= 0) 
        return (180/M_PI) * atan(_y / _x);//反正切，根据x,y所在象限确定角度具体计算公式

    if ( _x < 0 & _y >= 0) 
        return 180 - ( (180/M_PI) * atan(_y / (-_x)) );

    if ( _x < 0 & _y < 0) 
        return 180 + ( (180/M_PI) * atan(_y / _x) );

    if ( _x >= 0 & _y < 0)
        return 360 - ( (180/M_PI) * atan((-_y) / _x) );
} // xy2theta


MatrixXd circshift( MatrixXd &_mat, int _num_shift )//循环移位矩阵的列
{
    // shift columns to right direction 
    assert(_num_shift >= 0);

    if( _num_shift == 0 )
    {
        MatrixXd shifted_mat( _mat );
        return shifted_mat; // Early return 
    }

    MatrixXd shifted_mat = MatrixXd::Zero( _mat.rows(), _mat.cols() );//eigen库中的一个类，访问类的静态成员变量或函数，或命名空间中的成员。zero函数，参数为矩阵的行数和列数
    for ( int col_idx = 0; col_idx < _mat.cols(); col_idx++ )//遍历输入矩阵的每一列
    {
        int new_location = (col_idx + _num_shift) % _mat.cols();//取模。这行代码计算了新的列位置 new_location，以便在循环移位时将原始矩阵的列复制到正确的位置。
        shifted_mat.col(new_location) = _mat.col(col_idx);//shifted_mat 对象的 col() 方法，用于访问矩阵的某一列。将原始矩阵 _mat 的第 col_idx 列复制到新矩阵 shifted_mat 的新位置 new_location
    }

    return shifted_mat;

} // circshift


std::vector<float> eig2stdvec( MatrixXd _eigmat )//vec 的元素将与 _eigmat 中的元素相同，保持相同的顺序。目的是将 Eigen 矩阵的数据转换为一个浮点数向量
{
    std::vector<float> vec( _eigmat.data(), _eigmat.data() + _eigmat.size() );//浮点数数组容器。初始化。指向内部存储数组的指针。中元素的数量。
    return vec;
} // eig2stdvec


double SCManager::distDirectSC ( MatrixXd &_sc1, MatrixXd &_sc2 )
{
    int num_eff_cols = 0; // i.e., to exclude all-nonzero sector//即排除所有非0扇区
    double sum_sector_similarity = 0;
    for ( int col_idx = 0; col_idx < _sc1.cols(); col_idx++ )
    {
        VectorXd col_sc1 = _sc1.col(col_idx);//一种向量类型，一个动态大小的列向量。一维数组。.col(col_idx) 是 MatrixXd 类的成员函数，用于获取矩阵的第 col_idx 列作为一个 VectorXd。
        VectorXd col_sc2 = _sc2.col(col_idx);
        
        if( col_sc1.norm() == 0 | col_sc2.norm() == 0 )//某一列的向量范数为0
            continue; // don't count this sector pair. 

        double sector_similarity = col_sc1.dot(col_sc2) / (col_sc1.norm() * col_sc2.norm());//计算了两个向量的点积（内积）：col_sc1.dot(col_sc2)。将两个向量的对应元素相乘并求和。算了两个向量的范数（模）的乘积：col_sc1.norm() * col_sc2.norm()。范数是向量的长度或大小，通常使用L2 范数来计算。最后，我们将点积除以范数的乘积，得到一个相似度值。
//这行代码计算了两个向量之间的余弦相似度
        sum_sector_similarity = sum_sector_similarity + sector_similarity;
        num_eff_cols = num_eff_cols + 1;
    }
    
    double sc_sim = sum_sector_similarity / num_eff_cols;//所有有效列的平均相似度
    return 1.0 - sc_sim;

} // distDirectSC


int SCManager::fastAlignUsingVkey( MatrixXd & _vkey1, MatrixXd & _vkey2)//计算两个矩阵之间的对齐偏移量。找到使两个矩阵之间差异最小的偏移量
{
    int argmin_vkey_shift = 0;
    double min_veky_diff_norm = 10000000;
    for ( int shift_idx = 0; shift_idx < _vkey1.cols(); shift_idx++ )//对于每个偏移量执行操作
    {
        MatrixXd vkey2_shifted = circshift(_vkey2, shift_idx);//将 向右循环移位 列

        MatrixXd vkey_diff = _vkey1 - vkey2_shifted;//相减得到差异矩阵

        double cur_diff_norm = vkey_diff.norm();//差异矩阵的范数
        if( cur_diff_norm < min_veky_diff_norm )
        {
            argmin_vkey_shift = shift_idx;
            min_veky_diff_norm = cur_diff_norm;//更新
        }
    }

    return argmin_vkey_shift;

} // fastAlignUsingVkey


std::pair<double, int> SCManager::distanceBtnScanContext( MatrixXd &_sc1, MatrixXd &_sc2 )//返回类型，最小的扫描上下文距离和相应的移位索引
{
    // 1. fast align using variant key (not in original IROS18)使用可变键快速对齐(原始IROS18中没有)
    MatrixXd vkey_sc1 = makeSectorkeyFromScancontext( _sc1 );//函数计算 _sc1 和 _sc2 的变体键（variant key）。这些变体键用于快速对齐两个扫描上下文。
    MatrixXd vkey_sc2 = makeSectorkeyFromScancontext( _sc2 );
    int argmin_vkey_shift = fastAlignUsingVkey( vkey_sc1, vkey_sc2 );//函数找到一个初始的移位索引，以便在变体键上对齐 _sc1 和 _sc2。

    const int SEARCH_RADIUS = round( 0.5 * SEARCH_RATIO * _sc1.cols() ); // a half of search range .ratio头文件中设置为了0.1.  radius半径
    std::vector<int> shift_idx_search_space { argmin_vkey_shift };//初始化列表，用于在创建向量时，将 argmin_vkey_shift 的值作为初始元素添加到向量中。
    for ( int ii = 1; ii < SEARCH_RADIUS + 1; ii++ )
    {
        shift_idx_search_space.push_back( (argmin_vkey_shift + ii + _sc1.cols()) % _sc1.cols() );//这个计算式是将初始移位索引 argmin_vkey_shift 加上当前循环步数 ii，再加上 _sc1 矩阵的列数，然后对 _sc1 矩阵的列数取模。这表示向正方向（顺时针方向）进行移位
        shift_idx_search_space.push_back( (argmin_vkey_shift - ii + _sc1.cols()) % _sc1.cols() );//向负方向（逆时针方向）进行移位
    }
    std::sort(shift_idx_search_space.begin(), shift_idx_search_space.end());//排序函数。迭代器，指向向量的第一个元素
//std::sort 函数对 shift_idx_search_space 向量中的元素进行升序排序。
    // 2. fast columnwise diff 快速列式差分
    int argmin_shift = 0;//存储最小扫描上下文距离对应的移位索引。
    double min_sc_dist = 10000000;//存储最小的扫描上下文距离。
    for ( int num_shift: shift_idx_search_space )
    {
        MatrixXd sc2_shifted = circshift(_sc2, num_shift);//将 _sc2 矩阵循环移位 num_shift 个位置，并将结果存储在 sc2_shifted 中
        double cur_sc_dist = distDirectSC( _sc1, sc2_shifted );//计算扫描上下文之间的距离。
        if( cur_sc_dist < min_sc_dist )
        {
            argmin_shift = num_shift;
            min_sc_dist = cur_sc_dist;
        }
    }

    return make_pair(min_sc_dist, argmin_shift);
//括号 () 是用来调用函数的。make_pair 是一个函数，用于创建一个 std::pair 对象。参数。函数会返回一个包含这两个值的 std::pair 对象。
} // distanceBtnScanContext


MatrixXd SCManager::makeScancontext( pcl::PointCloud<SCPointType> & _scan_down )//从下采样的点云 _scan_down 创建一个扫描上下文
{
    TicToc t_making_desc;//计时

    int num_pts_scan_down = _scan_down.points.size();//获取下采样的点云 _scan_down 中的点的数量，并将其存储在 num_pts_scan_down 中。

    // main
    const int NO_POINT = -1000;//被用作一个标记，表示没有点
    MatrixXd desc = NO_POINT * MatrixXd::Ones(PC_NUM_RING, PC_NUM_SECTOR);//创建了一个 PC_NUM_RING 行 PC_NUM_SECTOR 列的矩阵 desc，并将所有元素初始化为 NO_POINT。传递给 Ones 函数的参数，表示创建的矩阵的行数和列数
//Eigen 库中的一个类，表示一个动态大小的矩阵。 C++ 中的作用域解析运算符，用于指定调用哪个类或命名空间的成员。Ones 是 MatrixXd 类的一个静态成员函数，用于创建一个元素全为 1 的矩阵。
    SCPointType pt;//被用来存储一个点。
    float azim_angle, azim_range; // wihtin 2d plane存储方位角和范围
    int ring_idx, sctor_idx;//存储环的索引和扇区的索引。
    for (int pt_idx = 0; pt_idx < num_pts_scan_down; pt_idx++)//pt表示point。遍历下采样的点云 _scan_down 中的所有点。
    {
        pt.x = _scan_down.points[pt_idx].x; 
        pt.y = _scan_down.points[pt_idx].y;
        pt.z = _scan_down.points[pt_idx].z + LIDAR_HEIGHT; // 获取当前点的 z 坐标，并加上激光雷达的高度。naive adding is ok (all points should be > 0).

        // xyz to ring, sector
        azim_range = sqrt(pt.x * pt.x + pt.y * pt.y);//计算当前点的方位范围
        azim_angle = xy2theta(pt.x, pt.y);//当前点的方位角

        // if range is out of roi, pass
        if( azim_range > PC_MAX_RADIUS )//当前点的方位范围是否超出了预设的最大半径。如果超出，就跳过当前点，处理下一个点。
            continue;

        ring_idx = std::max( std::min( PC_NUM_RING, int(ceil( (azim_range / PC_MAX_RADIUS) * PC_NUM_RING )) ), 1 );//计算当前点所在的环的索引。首先将当前点的方位范围 azim_range 除以预设的最大半径 PC_MAX_RADIUS，然后乘以环的数量 PC_NUM_RING。这样可以得到一个浮点数，表示当前点应该在哪个环上。将上面得到的浮点数向上取整，并转换为整数。这样可以确保环的索引是一个整数。
 //确保环的索引不会超过 PC_NUM_RING。如果上面的计算结果大于 PC_NUM_RING，则环的索引就是 PC_NUM_RING。确保环的索引不会小于 1。如果上面的计算结果小于 1，则环的索引就是 1。计算当前点所在的环的索引，索引的范围是从 1 到 PC_NUM_RING
        sctor_idx = std::max( std::min( PC_NUM_SECTOR, int(ceil( (azim_angle / 360.0) * PC_NUM_SECTOR )) ), 1 );//计算当前点所在的扇区的索引。
//int(ceil(2.3)) 的结果是 3，int(ceil(3.7)) 的结果是 4. ceil 是一个函数，用于向上取整.   int() 是一个类型转换操作，用于将一个值转换为整数。
        // taking maximum z 
        if ( desc(ring_idx-1, sctor_idx-1) < pt.z ) // -1 means cpp starts from 0访问描述符 desc 中特定位置的值。ring_idx-1 和 sctor_idx-1 是索引，分别表示环的索引和扇区的索引。在 C++ 中，索引是从 0 开始的，所以这里需要减 1。
            desc(ring_idx-1, sctor_idx-1) = pt.z; // update for taking maximum value at that bin当前点的 z 坐标 
    }//检查描述符 desc 中的值，如果该值小于当前点的 z 坐标 pt.z，则更新该值为 pt.z。

    // reset no points to zero (for cosine dist later)
    for ( int row_idx = 0; row_idx < desc.rows(); row_idx++ )// for 循环，用于遍历描述符 desc 的所有行。
        for ( int col_idx = 0; col_idx < desc.cols(); col_idx++ )//另一个 for 循环，嵌套在上一个循环中，用于遍历描述符 desc 的所有列
            if( desc(row_idx, col_idx) == NO_POINT )//括号 () 用于访问矩阵 desc 的元素。对象。索引。
                desc(row_idx, col_idx) = 0;

    t_making_desc.toc("PolarContext making");//记录了从 TicToc t_making_desc; 到现在的时间，并打印出来，标签为 “PolarContext making”

    return desc;//返回处理后的描述符 desc。
} // SCManager::makeScancontext


MatrixXd SCManager::makeRingkeyFromScancontext( Eigen::MatrixXd &_desc )//从扫描上下文 _desc 中创建一个环形键。环形键的每一行是 _desc 对应行的平均值。
{
    /* 
     * summary: rowwise mean vector
    */
    Eigen::MatrixXd invariant_key(_desc.rows(), 1);//行数为 _desc 的行数，列数为 1
    for ( int row_idx = 0; row_idx < _desc.rows(); row_idx++ )//遍历 _desc 的所有行。
    {
        Eigen::MatrixXd curr_row = _desc.row(row_idx);//获取 _desc 中的当前行 curr_row。
        invariant_key(row_idx, 0) = curr_row.mean();//计算 curr_row 的平均值，并将结果存储在 invariant_key 的相应位置。.mean() 是 Eigen 库中的一个成员函数，用于计算矩阵或向量的平均值，将 curr_row 中的所有元素相加，然后除以元素的数量。
    }

    return invariant_key;//返回环形键 
} // SCManager::makeRingkeyFromScancontext


MatrixXd SCManager::makeSectorkeyFromScancontext( Eigen::MatrixXd &_desc )//从扫描上下文 _desc 中创建一个扇区键。
{
    /* 
     * summary: columnwise mean vector
    */
    Eigen::MatrixXd variant_key(1, _desc.cols());//行数为 1，列数为 _desc 的列数。
    for ( int col_idx = 0; col_idx < _desc.cols(); col_idx++ )//遍历 _desc 的所有列。
    {
        Eigen::MatrixXd curr_col = _desc.col(col_idx);//获取 _desc 中的当前列 curr_col。
        variant_key(0, col_idx) = curr_col.mean();
    }

    return variant_key;//返回扇区键
} // SCManager::makeSectorkeyFromScancontext


void SCManager::makeAndSaveScancontextAndKeys( pcl::PointCloud<SCPointType> & _scan_down )
{
    Eigen::MatrixXd sc = makeScancontext(_scan_down); // v1 调用 makeScancontext 方法，从下采样的点云 _scan_down 中创建一个扫描上下文 sc。
    Eigen::MatrixXd ringkey = makeRingkeyFromScancontext( sc );//创建一个环形键 ringkey
    Eigen::MatrixXd sectorkey = makeSectorkeyFromScancontext( sc );//创建一个扇区键 sectorkey。
    std::vector<float> polarcontext_invkey_vec = eig2stdvec( ringkey );//将环形键 ringkey 转换为一个标准向量 polarcontext_invkey_vec

    polarcontexts_.push_back( sc ); //将扫描上下文 sc 添加到 polarcontexts_ 向量的末尾
    polarcontext_invkeys_.push_back( ringkey );
    polarcontext_vkeys_.push_back( sectorkey );
    polarcontext_invkeys_mat_.push_back( polarcontext_invkey_vec );

    // cout <<polarcontext_vkeys_.size() << endl;

} // SCManager::makeAndSaveScancontextAndKeys


std::pair<int, float> SCManager::detectLoopClosureID ( void )//检测循环闭合。循环闭合是指在构建点云地图时发现当前观测帧与之前某一帧具有相似的特征，并判断它们可能是同一个位置的不同观测。
{//类名，用于管理和处理点云地图中的各种操作。函数名，用于检测循环闭合。类中的成员函数。参数包括当前观测的极坐标描述符、最佳匹配帧索引和角度差。
    int loop_id { -1 }; // init with -1, -1 means no loop (== LeGO-LOAM's variable "closestHistoryFrameID")

    auto curr_key = polarcontext_invkeys_mat_.back(); // 将最后一个元素赋值给变量。当前观测帧的极坐标描述符对应的索引。current observation (query)
    auto curr_desc = polarcontexts_.back(); // current observation (query)//当前观测帧的极坐标描述符

    /* 
     * step 1: candidates from ringkey tree_
     */
    if( polarcontext_invkeys_mat_.size() < NUM_EXCLUDE_RECENT + 1)//元素数量
    {
        std::pair<int, float> result {loop_id, 0.0};//返回一个包含默认值(loop_id=0, dist=0.0) 的pair对象作为结果，并提前结束函数运行。
        return result; // Early return 
    }

    // tree_ reconstruction (not mandatory to make everytime)
    if( tree_making_period_conter % TREE_MAKING_PERIOD_ == 0) //取模。 to save computation cost
    {
        TicToc t_tree_construction;

        polarcontext_invkeys_to_search_.clear();//候选描述符。清空容器
        polarcontext_invkeys_to_search_.assign( polarcontext_invkeys_mat_.begin(), polarcontext_invkeys_mat_.end() - NUM_EXCLUDE_RECENT ) ;//赋值

        polarcontext_tree_.reset(); //重置
        polarcontext_tree_ = std::make_unique<InvKeyTree>(PC_NUM_RING /* dim */, polarcontext_invkeys_to_search_, 10 /* max leaf */ );//使用函数创建对象，调用对象的构造函数，传入三个参数。维度，待插入数据，最大叶子节点数量
        // tree_ptr_->index->buildIndex(); // inernally called in the constructor of InvKeyTree (for detail, refer the nanoflann and KDtreeVectorOfVectorsAdaptor)
        t_tree_construction.toc("Tree construction");//toc函数，计算  花费的时间
    }
    tree_making_period_conter = tree_making_period_conter + 1;
        
    double min_dist = 10000000; // init with somthing large可能被用来存储最小距离值，初始化为一个较大的数，以便在后续的计算中找到更小的值。
    int nn_align = 0;//通常，它可能被用来存储某种对齐或匹配的结果。
    int nn_idx = 0;//被用来存储最近邻（nearest neighbor）的索引

    // knn search//NUM_CANDIDATES_FROM_TREE 是一个整数，表示向量的大小。这个括号内的内容是向量的初始化参数，它告诉编译器我们想要创建一个多大的向量
    std::vector<size_t> candidate_indexes( NUM_CANDIDATES_FROM_TREE ); //初始化为括号内内容大小。定义了一个大小为 NUM_CANDIDATES_FROM_TREE 的 size_t 类型向量 candidate_indexes。这个向量可能被用来存储最近邻搜索的结果索引。
    std::vector<float> out_dists_sqr( NUM_CANDIDATES_FROM_TREE );//存储最近邻搜索的结果距离的平方。

    TicToc t_tree_search;//计时
    nanoflann::KNNResultSet<float> knnsearch_result( NUM_CANDIDATES_FROM_TREE );//定义了  类型的对象，并将其初始化为（）内内容。用于存储最近邻搜索的结果。
    knnsearch_result.init( &candidate_indexes[0], &out_dists_sqr[0] );//初始化，设置了存储搜索结果的位置和距离的向量。
    polarcontext_tree_->index->findNeighbors( knnsearch_result, &curr_key[0] /* query */, nanoflann::SearchParams(10) ); //这行代码在 polarcontext_tree_ 的索引中进行最近邻搜索，搜索的关键字是 curr_key，搜索参数是 nanoflann::SearchParams(10)。
//是使用 nanoflann 库在 polarcontext_tree_的索引中进行最近邻搜索。polarcontext_tree_->index：这是 polarcontext_tree_ 的索引，是一个KD树结构，用于存储和搜索数据。在KD树中进行最近邻搜索的函数，接受三个参数。是一个 nanoflann::KNNResultSet 对象，用于存储搜索结果。是搜索的关键字，也就是我们要找的点。是搜索参数，10 表示搜索时考虑的叶节点最大数量。
    t_tree_search.toc("Tree search");//这行代码记录了从 TicToc t_tree_search到现在的时间，并打印出来，标签为 “Tree search”。
//. 是 C++ 中的成员访问运算符，用于访问类、结构或联合的成员。toc 是 TicToc 类的一个成员函数，用于停止计时并打印出计时结果。括号 () 用于调用函数。在这里，toc 函数被调用，并且 "Tree search" 被传递给这个函数作为参数。调用 t_tree_search 对象的 toc 函数，并传入 "Tree search" 作为参数。这将停止 t_tree_search 的计时，并打印出标签为 “Tree search” 的计时结果。
    /* 
     *  step 2: pairwise distance (find optimal columnwise best-fit using cosine distance)
     */
    TicToc t_calc_dist;   
    for ( int candidate_iter_idx = 0; candidate_iter_idx < NUM_CANDIDATES_FROM_TREE; candidate_iter_idx++ )//for 循环，用于遍历所有的候选描述符。候选描述符的数量
    {//从 polarcontexts_ 数组中获取一个元素，并将其赋值给 polarcontext_candidate。存储了多个 MatrixXd 对象的数组。一个索引值，它从 candidate_indexes 数组中获取。当前的迭代次数
        MatrixXd polarcontext_candidate = polarcontexts_[ candidate_indexes[candidate_iter_idx] ];//获取了当前迭代的候选描述符 。访问向量 candidate_indexes 中索引为 candidate_iter_idx 的元素。
        std::pair<double, int> sc_dist_result = distanceBtnScanContext( curr_desc, polarcontext_candidate ); //调用了 distanceBtnScanContext 函数，计算了当前描述符 curr_desc 与候选描述符 polarcontext_candidate 之间的距离和对齐方式，结果存储在 sc_dist_result 中。
        
        double candidate_dist = sc_dist_result.first;//获取了 sc_dist_result 中的距离值 candidate_dist
        int candidate_align = sc_dist_result.second;//对齐方式

        if( candidate_dist < min_dist )//检查当前候选描述符的距离 candidate_dist 是否小于当前最小距离 min_dist。
        {
            min_dist = candidate_dist;//更新最小距离
            nn_align = candidate_align;//更新最近邻的对齐方式

            nn_idx = candidate_indexes[candidate_iter_idx];
        }//更新最近邻的索引
    }
    t_calc_dist.toc("Distance calc");//记录了从 TicToc t_calc_dist; 到现在的时间，并打印出来，标签为 “Distance calc”。

    /* 
     * loop threshold check
     */
    if( min_dist < SC_DIST_THRES )//最小距离。设定的阈值
    {
        loop_id = nn_idx; //最近邻索引
    
        // std::cout.precision(3); 
        cout << "[Loop found] Nearest distance: " << min_dist << " btn " << polarcontexts_.size()-1 << " and " << nn_idx << "." << endl;//打印出找到的循环，最近的距离，以及这个距离是在哪两个点之间的。
        cout << "[Loop found] yaw diff: " << nn_align * PC_UNIT_SECTORANGLE << " deg." << endl;//找到的循环的偏航角度差。
    }
    else
    {
        std::cout.precision(3); //设置打印的精度为3
        cout << "[Not loop] Nearest distance: " << min_dist << " btn " << polarcontexts_.size()-1 << " and " << nn_idx << "." << endl;
        cout << "[Not loop] yaw diff: " << nn_align * PC_UNIT_SECTORANGLE << " deg." << endl;
    }

    // To do: return also nn_align (i.e., yaw diff)
    float yaw_diff_rad = deg2rad(nn_align * PC_UNIT_SECTORANGLE);//首先将 nn_align（最近邻的对齐方式）乘以 PC_UNIT_SECTORANGLE（每个扇区的角度），然后将结果从度转换为弧度，得到偏航角度差 yaw_diff_rad。
    std::pair<int, float> result {loop_id, yaw_diff_rad};//定义了一个 std::pair 对象 result，并将其初始化为 {loop_id, yaw_diff_rad}。loop_id 是最近邻的索引，yaw_diff_rad 是偏航角度差

    return result;

} // SCManager::detectLoopClosureID

// } // namespace SC2