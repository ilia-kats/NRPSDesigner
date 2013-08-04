-- phpMyAdmin SQL Dump
-- version 3.4.11.1deb1
-- http://www.phpmyadmin.net
--
-- Server version: 5.5.31
-- PHP Version: 5.4.6-1ubuntu1.2

SET SQL_MODE="NO_AUTO_VALUE_ON_ZERO";
SET time_zone = "+00:00";

--
-- Database: `nrps_designer`
--
DROP DATABASE IF EXISTS `nrps_designer`;
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
  `origin_id` int(20) NOT NULL,
  `type_id` int(20) UNIQUE,
  `refseq_id` varchar(20) COLLATE utf8_unicode_ci DEFAULT NULL,
  `substrate_specificity_aa_id` int(20) DEFAULT NULL,
  `chirality` tinyint(1) NOT NULL,
  `uniprot_id` varchar(20) COLLATE utf8_unicode_ci default NULL,
  `description` text COLLATE utf8_unicode_ci NOT NULL,
  `pfam_linker_start` int(5) NOT NULL,
  `pfam_linker_stop` int(5) NOT NULL,
  `defined_linker_start` int(5) DEFAULT NULL,
  `defined_linker_stop` int(5) DEFAULT NULL,
  `pfam_start` int(5) NOT NULL,
  `pfam_stop` int(5) NOT NULL,
  `defined_start` int(5) DEFAULT NULL,
  `defined_stop` int(5) DEFAULT NULL,
  PRIMARY KEY (`domain_id`),
  KEY `substrate_specificity_aa_id` (`substrate_specificity_aa_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1;

--
-- Table structure for table `origin`
--

DROP TABLE IF EXISTS `origin`;
CREATE TABLE IF NOT EXISTS `origin` (
  `origin_id` int(20) NOT NULL AUTO_INCREMENT,
  `source` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `product` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `linkout` varchar(300) COLLATE utf8_unicode_ci NOT NULL,
  `dna_sequence` text COLLATE utf8_unicode_ci NOT NULL,
  `norine_id` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `description` text COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`origin_id`)
) ENGINE=InnoDB  DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;

-- --------------------------------------------------------

--
-- Table structure for table `substrates`
--

DROP TABLE IF EXISTS `substrates`;
CREATE TABLE IF NOT EXISTS `substrates` (
  `aa_id` int(20) NOT NULL UNIQUE AUTO_INCREMENT,
  `name` varchar(20) COLLATE utf8_unicode_ci NOT NULL,
  `structure` text COLLATE utf8_unicode_ci NOT NULL,
  `linkout` text COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`aa_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;

--
-- Table structure for table `types`
--

DROP TABLE IF EXISTS `types`;
CREATE TABLE IF NOT EXISTS `types` (
  `type_id` int(20) NOT NULL AUTO_INCREMENT,
  `name` varchar(4) COLLATE utf8_unicode_ci NOT NULL,
  `info` varchar(100) COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`type_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci AUTO_INCREMENT=1 ;


--
-- Table structure for table `origin_relations`
--
DROP TABLE IF EXISTS `origin_relations`;
CREATE TABLE IF NOT EXISTS `origin_relations` (
  `origin_basis` int(20) NOT NULL,
  `origin_derived` int(20) NOT NULL,
  `info` varchar(100) COLLATE utf8_unicode_ci NOT NULL,
  PRIMARY KEY (`origin_basis`, `origin_derived`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;


--
-- Constraints for tables `domain_types`
-- reference to main table
ALTER TABLE `domains`
  ADD CONSTRAINT `domains_ibfk_1` FOREIGN KEY (`origin_id`) REFERENCES `origin` (`origin_id`) ON DELETE CASCADE ON UPDATE CASCADE;
  
ALTER TABLE `domains`
  ADD INDEX ( `substrate_specificity_aa_id` ),
  ADD FOREIGN KEY (`substrate_specificity_aa_id`) REFERENCES `substrates` (`aa_id`) ON DELETE SET NULL ON UPDATE CASCADE;

ALTER TABLE `domains`
  ADD FOREIGN KEY (`type_id`) REFERENCES `types` (`type_id`) ON DELETE SET NULL ON UPDATE CASCADE;

ALTER TABLE `origin_relations`
  ADD FOREIGN KEY ( `origin_basis` ) REFERENCES `nrps_designer`.`origin` (`origin_id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD FOREIGN KEY ( `origin_derived` ) REFERENCES `nrps_designer`.`origin` (`origin_id`) ON DELETE CASCADE ON UPDATE CASCADE;