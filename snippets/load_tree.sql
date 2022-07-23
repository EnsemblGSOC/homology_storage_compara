SET @root = (SELECT root_id FROM gene_tree_root WHERE stable_id = 'ENSGT00730000111827');
SELECT * FROM gene_tree_node AS T1
LEFT JOIN (
	SELECT node_id, node_type, duplication_confidence_score FROM gene_tree_node_attr
) AS T2 ON T1.node_id = T2.node_id
LEFT JOIN (
	SELECT gene_member_id, seq_member_id FROM seq_member
) AS T3 ON T1.seq_member_id = T3.seq_member_id
LEFT JOIN (
	SELECT stable_id, gene_member_id FROM gene_member
) AS T4 ON T3.gene_member_id = T4.gene_member_id
WHERE root_id = @root;