SET @root = (SELECT root_id FROM gene_tree_root WHERE stable_id = 'ENSGT00730000111827');

SELECT *
FROM (
    SELECT * FROM
    (
        SELECT homology_id, gene_member_id AS gene_member_id_1
        FROM homology_member
    ) AS t1
    RIGHT JOIN(
        SELECT homology_id, gene_tree_node_id, description
        FROM homology
        WHERE gene_tree_root_id = @root
    ) AS t2
    USING (homology_id)
    LEFT JOIN (
        SELECT gene_member_id AS gene_member_id_1, stable_id AS stable_id_1
        FROM gene_member
    ) AS t3
    USING (gene_member_id_1)
) AS H1 JOIN
(
    SELECT * FROM
    (
        SELECT homology_id, gene_member_id AS gene_member_id_2
        FROM homology_member
    ) AS t4
    RIGHT JOIN(
        SELECT homology_id, gene_tree_node_id, description
        FROM homology
        WHERE gene_tree_root_id = @root
    ) AS t5
    USING (homology_id)
    LEFT JOIN (
        SELECT gene_member_id AS gene_member_id_2, stable_id AS stable_id_2
        FROM gene_member
    ) AS t6
    USING (gene_member_id_2)
) AS H2
USING (homology_id, gene_tree_node_id, description)
WHERE H1.gene_member_id_1 <> H2.gene_member_id_2;