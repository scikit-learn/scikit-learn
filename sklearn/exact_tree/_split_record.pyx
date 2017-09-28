cdef void split_record_reset(SplitRecord* split_record, int feature,
                                    int pos, float threshold,
                                    float impurity,
                                    float impurity_improvement, int nid,
                                    StatsNode* c_stats, StatsNode* l_stats,
                                    StatsNode* r_stats):
    split_record[0].feature = feature
    split_record[0].pos = pos
    split_record[0].threshold = threshold
    split_record[0].impurity = impurity
    split_record[0].impurity_improvement = impurity_improvement

    stats_node_copy_to(c_stats, &split_record.c_stats)
    stats_node_copy_to(l_stats, &split_record.l_stats)
    stats_node_copy_to(r_stats, &split_record.r_stats)


cdef inline void split_record_copy_to(SplitRecord* src_split_record,
                                      SplitRecord* dst_split_record):
    dst_split_record[0].feature = src_split_record[0].feature
    dst_split_record[0].pos = src_split_record[0].pos
    dst_split_record[0].threshold = src_split_record[0].threshold
    dst_split_record[0].impurity = src_split_record[0].impurity
    dst_split_record[0].impurity_improvement = src_split_record[0].impurity_improvement
    dst_split_record[0].nid = src_split_record[0].nid

    stats_node_copy_to(&src_split_record[0].c_stats,
                       &dst_split_record[0].c_stats)
    stats_node_copy_to(&src_split_record[0].l_stats,
                       &dst_split_record[0].l_stats)
    stats_node_copy_to(&src_split_record[0].r_stats,
                       &dst_split_record[0].r_stats)
