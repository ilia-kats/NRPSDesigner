-- phpMyAdmin SQL Dump
-- version 3.4.11.1deb1
-- http://www.phpmyadmin.net
--
-- Host: localhost
-- Generation Time: Jul 16, 2013 at 05:27 PM
-- Server version: 5.5.31
-- PHP Version: 5.4.6-1ubuntu1.2

SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";

--
-- Database: `nrps_designer`
--

CREATE DATABASE `nrps_designer` DEFAULT CHARACTER SET utf8 COLLATE utf8_unicode_ci;
USE `nrps_designer`;

-- --------------------------------------------------------
--
-- Table structure for table `domains`
--

DROP TABLE IF EXISTS `domains`;
CREATE TABLE IF NOT EXISTS `domains` (
  `domain_id` int(20) NOT NULL,
  `module_id` int(20) NOT NULL,
  `pathway_id` int(20) NOT NULL,
  `type_id` int(20) UNIQUE,
  `refseq_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `substrate_specificity_aa_id` int(20) default NULL,
  `chirality` tinyint(1) NOT NULL,
  `uniprot_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `bb_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `description` text COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_before` text COLLATE utf8_unicode_ci NOT NULL,
  `native_linker_after` text COLLATE utf8_unicode_ci NOT NULL,
  `pfam_border_before` int(5) NOT NULL,
  `pfam_border_after` int(5) NOT NULL,
  `defined_border_before` int(5) default NULL,
  `defined_border_after` int(5) default NULL,
  PRIMARY KEY (`domain_id`),
  KEY `substrate_specificity_aa_id` (`substrate_specificity_aa_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1;

--
-- Table structure for table `main`
--

DROP TABLE IF EXISTS `pathways`;
CREATE TABLE IF NOT EXISTS `pathways` (
  `pathway_id` int(20) NOT NULL AUTO_INCREMENT,
  `organism` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `pathway` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `linkout` varchar(300) COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text COLLATE utf8_unicode_ci NOT NULL,
  `uniprot_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `norine_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `description` text COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`pathway_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `substrates`
--

DROP TABLE IF EXISTS `substrates`;
CREATE TABLE IF NOT EXISTS `substrates` (
  `aa_id` int(20) NOT NULL UNIQUE AUTO_INCREMENT,
  `aminoacid` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `structure` text COLLATE utf8_unicode_ci NOT NULL,
  `mod_id` int(36),
  `linkout` text COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`aa_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;

--
-- Table structure for table `types`
--

DROP TABLE IF EXISTS `types`;
CREATE TABLE IF NOT EXISTS `types` (
  `type_id` int(20) NOT NULL AUTO_INCREMENT,
  `type_desc` varchar(100) COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`type_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;

--
-- Table structure for table `modifications`
--

DROP TABLE IF EXISTS `modifications`;
CREATE TABLE IF NOT EXISTS `modifications` (
  `modification_id` int(20) NOT NULL UNIQUE AUTO_INCREMENT,
  `modification_desc` varchar(100) COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`modification_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;

--
-- Table structure for table `substrates_xref_modifications`
--
DROP TABLE IF EXISTS `substrates_xref_modifications`;
CREATE TABLE IF NOT EXISTS `substrates_xref_modifications` (
  `aa_id` int(20) NOT NULL,
  `modification_id` int(20) NOT NULL,
  PRIMARY KEY (`aa_id`, `modification_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;


--
-- Constraints for tables `domain_types`
-- reference to main table
ALTER TABLE `domains`
  ADD CONSTRAINT `domains_ibfk_1` FOREIGN KEY (`pathway_id`) REFERENCES `pathways` (`pathway_id`) ON DELETE CASCADE ON UPDATE CASCADE;
  
ALTER TABLE `domains`
  ADD INDEX ( `substrate_specificity_aa_id` ),
  ADD FOREIGN KEY (`substrate_specificity_aa_id`) REFERENCES `substrates` (`aa_id`) ON DELETE SET NULL ON UPDATE CASCADE;

ALTER TABLE `domains`
  ADD FOREIGN KEY (`type_id`) REFERENCES `types` (`type_id`) ON DELETE SET NULL ON UPDATE CASCADE;

ALTER TABLE `substrates_xref_modifications`
  ADD FOREIGN KEY (`aa_id`) REFERENCES `substrates` (`aa_id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD FOREIGN KEY (`modification_id`) REFERENCES `modifications` (`modification_id`) ON DELETE CASCADE ON UPDATE CASCADE;