<html>
<head>
    <title>Insert NRPS</title>
    <style type="text/css">
        .valid {
            background-color: #FFFFFF;
        }
        .invalid {
            background-color: #FF0000;
        }
    </style>
</head>

<body>
    
<?php 

$host = "localhost"; // Host name 
$username = "root"; // Mysql username 
$password = ""; // Mysql password 
$db_name = "nrps_designer"; // Database name 
$domain_types_a_c = "domain_types_a_c"; // Table name 
echo $username;
// Connect to server and select database. 
mysql_connect("$host", "$username", "$password")or die("cannot connect"); 
mysql_select_db($db_name)or die("cannot select DB: ". mysql_error()); 



// Get values from form 
$pathway_id = mysql_escape_string($_POST['pathway_id']);
$genome_borders = mysql_escape_string($_POST['genome_borders']); 
$substrate_specificity = mysql_escape_string($_POST['substrate_specificity']);
$chirality = mysql_escape_string($_POST['chirality']);
$uni_prot_id = mysql_escape_string($_POST['uni_prot_id']);
$bb_id = mysql_escape_string($_POST['bb_id']);
$description = mysql_escape_string($_POST['description']);
$dna_sequence = mysql_escape_string($_POST['dna_sequence']);
$native_linker_before = mysql_escape_string($_POST['native_linker_before']);
$native_linker_after = mysql_escape_string($_POST['native_linker_after']);

$invalid = false;
if (!$pathway_id && !$genome_borders && !$substrate_specificity && !$chirality && !$uni_prot_id && !$bb_id) {
    $clear_form = true;
} else {
    if (!$pathway_id) {
        $pathway_id_invalid = true;
        $invalid = true;
    }
    if (!$genome_borders) {
        $genome_borders_invalid = true;
        $invalid = true;
    }
    if (!$substrate_specificity) {
        $substrate_specificity_invalid = true;
        $invalid = true;
    }
    if (!$chirality) {
        $chirality_invalid = true;
        $invalid = true;
    }
    if (!$uni_prot_id) {
        $uni_prot_id_invalid = true;
        $invalid = true;
    }
    if (!$bb_id) {
        $bb_id_invalid = true;
        $invalid = true;
    }
    if (!$description) {
        $description_invalid = true;
        $invalid = true;
    }
    if (!$dna_sequence) {
        $dna_sequence_invalid = true;
        $invalid = true;
    }
    if (!$native_linker_before) {
        $native_linker_before_invalid = true;
        $invalid = true;
    }
    if (!$native_linker_after) {
        $native_linker_after_invalid = true;
        $invalid = true;
    }
}

if (!$invalid) {
             // Insert data into mysql 
$sql="INSERT INTO $domain_types_a_c SET
`pathway_id` = '$pathway_id',
`genome_borders` = '$genome_borders',
`substrate_specificity_aa_id` = '$substrate_specificity',
`chirality` = '$chirality',
`uni_prot_id` = '$uni_prot_id',
`bb_id` = '$bb_id',
`description` ='$description',
`dna_sequence`= '$dna_sequence',
`native_linker_before`= '$native_linker_before',
`native_linker_after`= '$native_linker_after'";
$result = mysql_query($sql); 

// if successfully insert data into database, displays message "Successful". 
if($result){ 
echo "Successful"; 
echo "<BR>"; 
echo "<a href='insert3.php'>Back to main page</a>";
$clear_form = true;
} else { 
    echo mysql_error(); 
} 
    
}
// close connection 
mysql_close();
function getClass($invalid) {
    return ($invalid ? "invalid" : "valid");
}
function getValue($var) {
    global $clear_form;
    return ($clear_form ? "" : $_POST[$var]);
}
?>
<form method="post" action="<?php echo $_SERVER["PHP_SELF"];?>">
        Genome_Borders: <input type="text" name="organism" class="<?php echo getClass($pathway_id_invalid);?>" value="<?php echo getValue('organism');?>"/><br/>
        Pathway_ID: <input type="text" name="pathway_id" class="<?php echo getClass($genome_borders_invalid);?>" value="<?php echo getValue('pathway');?>"/><br/>
        Substrate_Specificity: <input type="text" name="linkout" class="<?php echo getClass($substrate_specificity_invalid);?>" value="<?php echo getValue('substrate_specificity');?>"/><br/>
        dna_sequence: <input type="text" name="dna_sequence" class="<?php echo getClass($chirality_invalid);?>" value="<?php echo getValue('dna_sequence');?>"/><br/>
        UniProt_ID: <input type="text" name="uni_prot_id" class="<?php echo getClass($uni_prot_id_invalid);?>" value="<?php echo getValue('uni_prot_id');?>"/><br/>
        Biobrick_ID: <input type="text" name="bb_id" class="<?php echo getClass($bb_id_invalid);?>" value="<?php echo getValue('bb_id');?>"/><br/>
        Description: <input type="text" name="description" class="<?php echo getClass($description_invalid);?>" value="<?php echo getValue('description');?>" /><br/>
        DNA_Sequence: <input type="text" name="dna_sequence" class="<?php echo getClass($dna_sequence_invalid);?>" value="<?php echo getValue('dna_sequence');?>" /><br/>
        Native_Linker_Before: <input type="text" name="native_linker_before" class="<?php echo getClass($native_linker_before_invalid);?>" value="<?php echo getValue('native_linker_before');?>" /><br/>
        Native_Linker_After: <input type="text" name="native_linker_after" class="<?php echo getClass($native_linker_after_invalid);?>" value="<?php echo getValue('native_linker_after');?>" /><br/>
        <input type="submit" value="Insert into database"/>
</form>


</body>
</html>
