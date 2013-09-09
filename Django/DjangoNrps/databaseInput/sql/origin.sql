DROP FUNCTION IF EXISTS `get_origin_parent`;
CREATE DEFINER = CURRENT_USER FUNCTION `get_origin_parent`(`value` INT) RETURNS INT(11) NOT DETERMINISTIC READS SQL DATA SQL SECURITY INVOKER
BEGIN
    DECLARE _id INT;
    DECLARE CONTINUE HANDLER FOR NOT FOUND SET @id = NULL;

    IF @id IS NULL THEN
            RETURN NULL;
    END IF;
    SET _id = @id;

    SELECT  parent_id
    INTO    @id
    FROM    databaseInput_origin
    WHERE   id = @id;
    SET @level := @level + 1;
    RETURN _id;
END

DROP PROCEDURE IF EXISTS `get_origin_hierarchy`;
CREATE DEFINER = CURRENT_USER PROCEDURE `get_origin_hierarchy`(IN `id` INT) NOT DETERMINISTIC READS SQL DATA SQL SECURITY INVOKER
SELECT  ori.*
FROM    (
    SELECT  get_origin_parent(id) AS id, @level AS level
    FROM    (
        SELECT  @id := id, @level := 0
    ) vars, databaseInput_origin
) vars
JOIN    databaseInput_origin ori
ON      ori.id = vars.id ORDER BY level DESC;